function [Q0, q0, Qs, qs, b, x0, f_fun, G_fun] = QCQP(m, n)
% This function generates the data (Q0, q0, Qs, qs, b, x0) for
% a convex quadratically constrained quadratic program (QCQP)

% Sparsity level
sp = .01;

% Create positive semidefinite Q0
U = sprand(n,n,sp);
d = 100*rand(n, 1);
Q0 = U'*diag(d)*U;
Q0 = .5*(Q0'+Q0) + 1e-13*eye(n);
q0 = 10+randn(n, 1);

% Create cell arrays of positive semidefinite 
% matrices Qs{i} and vectors qs{i}
Qs = cell(m, 1);
qs = cell(m, 1);
for i=1:m
    U = sprand(n,n,sp);
    d = 100*rand(n, 1);
    Qs{i} =  U'*diag(d)*U;
    Qs{i} = .5*(Qs{i} + Qs{i}') + 1e-13*eye(n);
    qs{i} = 10+randn(n, 1);
end

% x0=0 is a strictly feasible point
b = 20*ones(m, 1);
x0 = zeros(n, 1);

f_fun = @(x) QuadObj(x, Q0, q0);
G_fun = @(x) QuadConst(x, Qs, qs, b);
end

% Define quadratic constraints
function [Gx, DGx] = QuadConst(x, Qs, qs, b)
m = length(b);
n = length(x);
Gx = zeros(m, 1);
if nargout>1
    DGx = zeros(m, n);
end
for i=1:m
    Qx = Qs{i}*x; 
    Gx(i) = x'*(1/2*Qx + qs{i}) - b(i);
    if nargout>1
        DGx(i, :) = Qx + qs{i};
    end
end
end

% Define the quadratic part of objective
function [fx, dfx] = QuadObj(x, Q, q)
Qx = Q*x; 
fx = x'*(1/2*Qx + q);
if nargout>1
    dfx = Qx + q;
end
end