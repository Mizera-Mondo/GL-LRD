function [P, Q] = LRD_normed(X, k, options)
%LRD_NORMED Normalized Low-Rank Decomposition: X = P*Q' with max rank = k.
% ||X - PQ'||_F^2 + ita*(||P||_F^2 + ||Q||_F^2)
% ita is the normalization parameter.
% Developed to deal with non-convexity of LRD. (X = PQ' is basically NOT
% unique, ref to SVD.)
arguments
    X double
    k double
    options.ita = 1
    options.solver = 'CVX'
end
ita = options.ita;
[a, b] = size(X);
P = randn(a, k);
Q = randn(b, k);

if strcmp(options.solver, 'CVX')

elseif strcmp(options.solver, 'quadprog')

else
    error('%s is not a vaild solver.', options.solver);
end
end

