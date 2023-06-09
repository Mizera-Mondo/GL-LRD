function [X, L] = GL_NO_LRD(Y, R, k, alpha, beta)
%LR-DGI Solve 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} +
%beta/2||L||_F^2
% k is the maximum rank of estimated signal matrix X

global debug;
debug = true;

X = Y;


% Random vaild initial value of A
[n, T] = size(Y);
B = [zeros(T - 1, 1) eye(T - 1); ...
    zeros(1, T)];
D = @(X) X - R*X*B;
A = rand(n, n);
A = A - diag(diag(A));
A = A + A';
A = n*A./sum(A, 'all');
L = diag(sum(A)) - A;

tol = 1e-2;
iter = 1;
maxIter = 1000;
isConverge = false;
isMaxIter = false;

if debug
    disp('=========================================================');
    disp('Starting Graph Learning. Low-rank decomposition disabled.');
    disp('=========================================================');

end

while ~isConverge && ~isMaxIter
    
    if debug
        disp(['Current Iteration: ' num2str(iter)]);
    end

    L_old = L;
    X_old = X;
    % Optimizing L a.k.a. A
    DX = D(X);
    M1 = DX*DX';
    M2 = repmat(diag(M1), 1, n);
    M = M2 + M2' - 2*M1;
    A = solveSubA(M, alpha, beta);

    L = diag(sum(A)) - A;

    % Optimizing X
    X = solveSubX(Y, L, R, B, alpha, k);

    isConverge = (norm(L_old - L, 'fro')/(norm(L_old, 'fro') + 1e-10) < tol) && (norm(X_old - X, 'fro')/(norm(X_old, 'fro') + 1e-10) < tol);
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
    
end
end

function A = solveSubA(M, alpha, beta)
%solveSubA solve the sub-problem ||A||_F^2 + alpha/beta*Tr{AM}
[n, ~] = size(M);
on = ones(n, 1);
ze = zeros(n, 1);
cvx_begin quiet
    variable A(n, n) symmetric nonnegative
    minimize square_pos(norm(A, 'fro')) + alpha/beta*trace(A*M)
    subject to
        diag(A) == ze;
        on'*A*on == n;
cvx_end
end

function X = solveSubX(Y, L, R, B, alpha, k)
%solveSubX solve the sub-problem 1/2||D(Y-X)||_F^2 +
%alpha*Tr{D(X)'*L*D(X)}, s.t. rank(X) <= k, using ADMM

[n, ~] = size(Y);
% P = randn(n, k); % random initial value
% Q = (P\Y)';
X = Y;
% Th = X - P*Q';
rho = 1;
ita = 1.05;
rhoMax = 10;

tolSqu = 1e-4;
tol = sqrt(tolSqu);
maxIter = 1000;
iter = 1;
isConverge = false;
isMaxIter = false;

while ~isConverge && ~isMaxIter
    X_old = X;
    % P_old = P;
    % Q_old = Q;

    % Update of X
    % X_ = P*Q' - Th/rho;
    X = updateX(X, Y, L, R, B, alpha, rho);

    % Update of P
    %P = (X + 1/rho*Th)*Q/(Q'*Q);
    % Update of Q
    %Q = (X' + 1/rho*Th')*P/(P'*P);
    % Update of Th and rho
    %Th = Th + rho*(X - P*Q');
    %rho = min([ita*rho, rhoMax]);

    % Terminating Condition Check
    if norm(X_old - X, 'fro')/(norm(X_old, 'fro') + 1e-10) < tol
        isConverge = true;
    end
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
end

function X = updateX(X, Y, L, R, B, alpha, rho)
%updateX solves the problem: 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} 

D = @(X) X - R*X*B;
DY = D(Y);
L_ = eye(size(L)) + 2*alpha*L;
H = DY - R'*DY*B';

tarFun = @(X) 1/2*norm(DY - D(X), 'fro')^2 + alpha*trace(D(X)'*L*D(X));
grad = @(X) L_*X + R'*L_*R*X*(B*B') - L_*R*X*B - R'*L_*X*B' + rho*X - H;

% TODO: Armijo update
c = 1e-2;
a = 0.9;
iter = 1;
maxIter = 100;
isArmijoNod = false;
isMaxIter = false;

gradX = grad(X);
while ~isArmijoNod && ~isMaxIter
    deltaX = -1*a.^(iter)*gradX;
    expectDescent = c*trace(gradX'*deltaX);
    realDescent = tarFun(X + deltaX) - tarFun(X);
    isArmijoNod = realDescent < expectDescent;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
if isArmijoNod
    X = X + deltaX;
else
%    disp('X unchanged due to non-decreasing within tolerance.');
end
end