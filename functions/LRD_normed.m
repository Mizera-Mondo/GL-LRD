function [P, Q] = LRD_normed(X, k, options)
%LRD_NORMED Normalized Low-Rank Decomposition: X = PQ' with max rank = k.
% ||X - PQ'||_F^2 + ita*(||P||_F^2 + ||Q||_F^2)
% ita is the normalization parameter.
% Developed to deal with non-convexity of LRD. (X = PQ' is basically NOT
% unique, ref to SVD.)
arguments
    X double
    k double
    options.ita = 1
    options.solver = 'CVX'
    options.maxIter = 1000
    options.tolerance = 1e-3
end
ita = options.ita;
tol = options.tolerance;
[a, b] = size(X);
P = randn(a, k);
Q = randn(b, k);
I = eye(k);
if strcmp(options.solver, 'CVX')
    iter = 1;
    isConverge = false;
    isMaxIter = false;
    while ~isConverge && ~isMaxIter
        P_old = P;
        Q_old = Q;

        cvx_begin quiet
        variables P(a, k)
        minimize square_pos(norm(X - P*Q', 'fro')) + ita*square_pos(norm(P, 'fro'))
        cvx_end

        cvx_begin quiet
        variables Q(b, k)
        minimize square_pos(norm(X - P*Q', 'fro')) + ita*square_pos(norm(Q, 'fro'))
        cvx_end

        isConverge = norm(P - P_old, 'fro')/norm(P_old) < tol && ...
            norm(Q - Q_old, 'fro')/norm(Q_old) < tol;
        isMaxIter = iter >= options.maxIter;
        iter = iter + 1;
    end


elseif strcmp(options.solver, 'closedform')
    % error('%s is under constrution and not available yet.', options.solver);
    iter = 1;
    isConverge = false;
    isMaxIter = false;
    while ~isConverge && ~isMaxIter
        P_old = P;
        Q_old = Q;

        P = X*Q/(ita*I + Q'*Q);
        Q = X'*P/(ita*I + P'*P);

        isConverge = norm(P - P_old, 'fro')/norm(P_old) < tol && ...
            norm(Q - Q_old, 'fro')/norm(Q_old) < tol;
        isMaxIter = iter >= options.maxIter;
        iter = iter + 1;
    end


else
    error('%s is not a vaild solver.', options.solver);
end
end

