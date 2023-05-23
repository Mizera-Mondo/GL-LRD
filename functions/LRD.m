function [P, Q] = LRD(X, k)
%LRD Low-rank decomposition of matrix X = P*Q' with maximum rank K with ADMM
[n, T] = size(X);
P = randn(n, k);
Q = randn(T, k);
Th = X - P*Q';

rho = 1;
ita = 1.05;
rhoMax = 10;
tol = 1e-3;

iter = 1;
maxIter = 1000;
isAdmmConverge = false;
isMaxIter = false;

while ~isAdmmConverge && ~isMaxIter
    P_old = P;
    Q_old = Q;

    cvx_begin quiet
        variable P(n, k)
        minimize trace(Th'*(X - P*Q')) + rho/2*square_pos(norm(X - P*Q', 'fro'))
    cvx_end

    cvx_begin quiet
        variable Q(T, k)
        minimize trace(Th'*(X - P*Q')) + rho/2*square_pos(norm(X - P*Q', 'fro'))
    cvx_end
    Th = Th + rho*(X - P*Q');
    rho = min([ita*rho, rhoMax]);
    isAdmmConverge = norm(P - P_old, 'fro')/norm(P_old) < tol ...
                  && norm(Q - Q_old, 'fro')/norm(Q_old) < tol;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;

end

