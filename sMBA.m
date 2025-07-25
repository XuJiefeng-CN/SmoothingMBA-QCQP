function [x, psi, g, out] = sMBA(Q0, q0, QQ, qq, b, rho1, K, opt)
% initialization
x = opt.x0;         % initial point  
rbar = opt.rbar;

ls = zeros(K+1,2);    % total iteration of line-search
psi_s = [];         % recording the values of objective f
gs = [];            % g = max_i g_i
Lf = 1;
Lg = 1;

% default parameters
alpha4 = .01;

tau1 = .01;   
tau2 = 1e-6; 
L_check = 1e-8;
L_hat = 1e8;

n0 = 300;
nu0 = 1/(10*n0 + 1);
r0 = .1;

% cell to array
n = length(x);
QQ = cell2mat(QQ); % size of n*m x n 
QQ = reshape(QQ, n, []); % size of n x m*n
qq = cell2mat(qq);
qq = reshape(qq, n, [])'; % size of m x n

% require g(x0)<0
Qx = x'*QQ;                 % size of 1 x m*n
Qx = reshape(Qx, [], n);    % size of m x n
G = (Qx/2 + qq)*x - b;
g = max(G);
if g>0
    disp('g>0')
    psi = [];
    out = [];
    return
end

% find mu0>0 s.t. g_mu0(x0)<.1*g(x0)
mu0 = .5;      % smoothing parameter     
while 1
    % g_mu0(x0)
    u = exp((G-g)/mu0);
    u_sum = sum(u);
    g_mu = g + mu0*log(u_sum) + alpha4*mu0;

    % break condition 
    if g_mu<= .1*g
        break
    end

    % decrease mu0
    mu0 = .5*mu0;
end
mu = mu0;

% psi(x0)
Q0x = Q0*x; 
f = x'*(Q0x/2 + q0);
df = Q0x + q0;
psi = f + rho1*sum(abs(x));

%  nabla g_mu0(x0)
dh_mu = u/u_sum;            % nabla h_mu(G(x))
DG = Qx + qq;
dg_mu = (dh_mu'*DG)';       % nabla g_mu(x)

psi_s = [psi_s; psi];
gs = [gs; g];

% print
fprintf('---------------------------------------------------------------------------------\n');
fprintf('%5s   %5s   %5s   %10s   %9s   %8s   %8s   %8s\n', ...
        'k', 'is', 'js', 'psi', 'g', 'mu', 'Lf', 'Lg');
fprintf('---------------------------------------------------------------------------------\n');

for k = 0:K
    % line-seach
    i = 0;
    for j = 0:40
        Lmu = Lg/mu;
        % Solve subproblem
        RR = (dg_mu/Lmu)'*(dg_mu/Lmu) - 2*g_mu/Lmu;
        xbar = x - 1/Lf * df;
        xhat = x - 1/Lmu * dg_mu;
        [x1, lambda1] = SubP_alpha(xbar, xhat, RR, Lf/rho1);

        % calculate Qx1, Gx1 and gx1
        Qx = x1'*QQ;                 % size of 1 x m*n
        Qx = reshape(Qx, [], n);    % size of m x n
        G = (Qx/2 + qq)*x1 - b;
        g = max(G);

        % calculate g1_mu := g_mu(x1)
        u = exp((G-g)/mu);
        u_sum = sum(u);
        g1_mu = g + mu*log(u_sum) + alpha4*mu;
        
        if g1_mu<=0      % is feasible
            % calculate psi(x1) = f(x1) + P1(x1)
            Q0x1 = Q0*x1; 
            f1 = x1'*(Q0x1/2 + q0);
            df1 = Q0x1 + q0;
            psi1 = f1 + rho1*sum(abs(x1));
            s = x1 - x;
            ss = s'*s;

            % Sufficient decreasing condition
            is_suf_dec = psi1 <= psi - .5*(tau1 + tau2*lambda1/mu)*ss;
            if is_suf_dec
                % Both conditions hold
                break
            else
                % increase Lf
                Lf = 2*Lf;
                i = i+1;
            end
        end

        % increase Lg
        Lg = 2*Lg;
    end
    ls(k+1, 1) = i;
    ls(k+1, 2) = j;

    if j>=40 
        disp('j>40, infeasible!!!')
        break
    end

    % update mu(k+1)
    r = r0 + (k+1)/K*(rbar-r0);     % r(k+1)
    k1 = mod(k+1, n0+1);            % remainder
    kprod = k+1 - k1;
    kbar = kprod + nu0*k1 + 1;
    mu = mu0*kbar^(-r);

    % calculate dg_mu1
    u = exp((G-g)/mu);
    u_sum = sum(u);
    g_mu = g + mu*log(u_sum) + alpha4*mu;
    dh_x1mu1 = u/u_sum;            % nabla h_mu
    DG = Qx + qq;
    dg_x1mu1 = (dh_x1mu1'*DG)';          % nabla g_mu

    % calculate BB stepsize: Lf_bb
    delta_f = df1 - df;
    Lf = max(L_check, Lf/2);
    if sqrt(ss)>1e-12
        Lf_BB = abs(delta_f'*s/ss);      % BB stepsize of f
        if Lf_BB >= L_check && Lf_BB <= L_hat
            Lf = Lf_BB;
        end
    end
    
    % calculate BB stepsize: Lg_bb
    delta_g = dg_x1mu1 - dg_mu;
    Lg = max(L_check, Lg/2);
    if sqrt(ss)>1e-12
        Lmu_bb = abs(delta_g'*s/ss);
        Lg_bb = Lmu_bb*mu;            % BB stepsize of gmu
        if Lg_bb >= L_check && Lg_bb <= L_hat
            Lg = Lg_bb;
        end
    end

    % update x, psi and g
    x = x1;
    df = df1;
    psi = psi1;
    dg_mu = dg_x1mu1;

    psi_s = [psi_s; psi];
    gs = [gs; g];
    
    % Print
    if mod(k, 500)==0
        fprintf('%5d | %5d | %5d | %10.5f | %9.2e | %8.2e | %8.2e | %8.2e\n', ...
                    k, sum(ls(:, 1)), sum(ls(:, 2)), psi, g, mu, Lf, Lg);
    end
end      
fprintf('---------------------------------------------------------------------------------\n\n');

% output
out.k = k;
out.mu0 = mu0;
out.mu = mu;
out.L_mu = Lmu;
out.ls = ls;
out.psi_vals = psi_s;
out.gvals = gs;

end