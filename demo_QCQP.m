%% Demo: Convex QCQP with ℓ₁-Regularization
% Problem formulation:
%   Objective:    ψ(x) := f(x) + ρ * ||x||₁, 
%                 where f(x) = x' * (0.5 * Q0 * x + q0)
%   Constraints:  gᵢ(x) := x' * (0.5 * Qs{i} * x + qs{i}) - bᵢ ≤ 0,  for i = 1,...,m
clear; clc

% Problem identifier
prob = 'a';

% Set up a results recorder
timestamp = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
foldername = 'results';
if ~exist(foldername, 'dir')
    mkdir(foldername);
end
fID = fopen([foldername '/results.txt'], 'a');

%% Switch problem
K = 5000;   % Maximum number of iterations
switch prob
    case 'a'
        % problem
        rand('seed', 2024)
        randn('seed', 2024)
        n = 500; m = 100;
        [Q0, q0, Qs, qs, b, x0, f_fun, G_fun] = QCQP(m, n);
        rho = 1;

        % Decreasing rate of mu_k
        rbars = [.51 0.75 .99]; 

    case 'b'
        % problem
        rand('seed', 2025)
        randn('seed', 2025)
        n = 1000; m = 500;
        [Q0, q0, Qs, qs, b, x0, f_fun, G_fun] = QCQP(m, n);
        rho = 1;

        % Decreasing rate of mu_k
        rbars = [.51 0.75 .99]; 
end

%% print
fprintf(fID, '  Ell-one regularized QCQP with (n,m) = (%3d,%3d), runing at %10s \n', n, m, char(timestamp));
fprintf(fID, '%5s   %8s   %9s   %8s   %5s   %7s   %4s   %8s   %8s\n', ...
             'rbar', 'f', 'g', 'time', 'iter', 'mean ls', 'mu0', 'mu', 'L_mu');
fprintf(fID, '-------------------------------------------------------------------------------------------\n');
fprintf('  Ell-one regularized QCQP with (n,m) = (%3d,%3d), runing at %10s \n', n, m, char(timestamp));
fprintf('-------------------------------------------------------------------------------------------\n');

%% cvx solver
fprintf('cvx...\n')
fprintf('-------------------------------------------------------------------------------------------\n');
st = tic;
cvx_solver SDPT3
cvx_begin
    variable x_cvx(n)
    minimize( .5* quad_form(x_cvx, Q0) + q0' * x_cvx + rho*norm(x_cvx, 1) )
    subject to
    for i=1:m
        .5* quad_form(x_cvx, Qs{i}) + qs{i}' * x_cvx <= b(i);
    end
cvx_end
time_cvx = toc(st);
g_cvx = max(G_fun(x_cvx));

% save results
fprintf(fID, '%5s | %10.5f | %9.2e | %8.5f\n', ...
    'cvx', cvx_optval, g_cvx, time_cvx);

%% sMBA methods
for t = 1:length(rbars)
    rbar = rbars(t);
    Mths{t} = sprintf('(%1d) $\bar r = %.2g$', t, rbar);
    fprintf('rbar = %5.2g...\n', rbar)

    % options
    opt = [];           
    opt.x0 = x0;
    opt.rbar = rbar; 
    
    % run sMBA
    tic_temp = tic;
    [x{t}, psi_val{t}, gval{t}, out{t}] = ...
        sMBA(Q0, q0, Qs, qs, b, rho, K, opt);
    time(t) = toc(tic_temp);
    
    % save results
    fprintf(fID, '%5.2g | %10.5f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e\n', ...
        rbars(t), psi_val{t}, gval{t}, time(t), ...
        out{t}.k, mean(out{t}.ls(:, 2),'all'), out{t}.mu0, out{t}.mu, out{t}.L_mu);
end
tt = t;
fprintf(fID, '-------------------------------------------------------------------------------------------\n');

%% print Latex table
style_tex = {'{\color{blue}\bf - - - }', '{\color{blue}\bf $\cdots$ }', '{\color{blue}\bf --- }'};
fprintf(fID, '                       tabale (%s)      \n', prob);
fprintf(fID, '----------------------------------------------------------------------\n');
fprintf(fID, 'line style  &  $\\bar r$ & objective & \\text{time (s)} \\\\\n');
for t=1:tt
%     line style & item &  $\bar r$& $\psi^{K}$ & \text{time (s)} \\
    fprintf(fID, '%27s & $%4.2g$ & $%12.5f$ & $%7.1f$\\\\\n',...
        style_tex{t}, rbars(t), psi_val{t}, time(t));
end
fprintf(fID, '%27s &  %4s  & $%12.5f$ & $%7.1f$ \n', '', 'cvx', cvx_optval, time_cvx);
fprintf(fID, '----------------------------------------------------------------------\n\n');
fclose(fID);

%% figure
figure; hold on
pos = [4 1.5 4 3];
set(gcf,'Units','Inches');
set(gcf,'Position',pos);
style = {'--b', ':b', '-b'};

% plot sMBA
for t = 1:tt
    temp = (out{t}.psi_vals-cvx_optval)/max(1,abs(cvx_optval));
    plot(log10(temp(2:end-1)), style{t});
    methods{t} = sprintf('rbar=%.3g', rbars(t));
end

% label
title(sprintf('(%s) $(n, m) = (%d, %d)$', prob, n, m), 'Interpreter','latex')
xlabel('iteration $k$','Interpreter','latex')
ytickformat('10^{%.2g}')
ylabel('$\log \omega_{k} $', 'Interpreter','latex')
legend(methods)

% figure to pdf
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, [foldername '/' prob '.pdf'], '-dpdf', '-r0');
