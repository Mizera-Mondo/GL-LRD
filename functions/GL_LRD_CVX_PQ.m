function [X, L] = GL_LRD_CVX_PQ(Y, R, k, alpha, beta)
%LR-DGI Solve 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} +
%beta/2||L||_F^2
% k is the maximum rank of estimated signal matrix X
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

iter = 1;
maxIter = 1000;
isConverge = false;
isMaxIter = false;

while ~isConverge && ~isMaxIter
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

    isConverge = norm(L_old - L, 'fro')/(norm(L_old, 'fro') + 1e-10) < tol ...
            && norm(X_old - X, 'fro')/(norm(X_old, 'fro') + 1e-10) < tol;
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

%% For Debugging
D = @(X) X - R*X*B;
tarFun = @(X) 1/2*norm(D(Y) - D(X), 'fro')^2 + alpha*trace(D(X)'*L*D(X));
cstrFun = @(X, P, Q) norm(X - P*Q', 'fro');
largFun = @(X, P, Q, Th, rho) tarFun(X) + rho/2*cstrFun(X, P, Q) + trace(Th'*(X - P*Q'));
debug = false;
%%
[n, T] = size(Y);
P = randn(n, k); % random initial value
Q = (P\Y)';
X = Y;
Th = X - P*Q';
rho = 0.05;
ita = 1.05;
rhoMax = 10;

tolSqu = 1e-6;
tol = sqrt(tolSqu);
maxIter = 1000;
iter = 1;
isConverge = false;
isMaxIter = false;

while ~isConverge && ~isMaxIter
    % Debug output
    if debug
        disp(['Iter: ' num2str(iter)]);
        disp(['Target Function starts at ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation starts at : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
    end
    X_old = X;
    P_old = P;
    Q_old = Q;
    % Update of X

    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: Updating X');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end   


    X_ = P*Q' - Th/rho;
    X = updateX(X, Y, L, X_, R, B, alpha, rho);

    % Update of P and Q with CVX
    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: Updating P');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end       
    cvx_begin quiet
        variable P(n, k)
        minimize rho/2*norm(X - P*Q' + Th/rho, 'fro')
    cvx_end

    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: Updating Q');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end   

    cvx_begin quiet
        variable Q(T, k)
        minimize rho/2*norm(X - P*Q' + Th/rho, 'fro')
    cvx_end

    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: Updating Th');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end 

    % % Update of P
    % P = (X + 1/rho*Th)*Q/(Q'*Q);
    % % Update of Q
    % Q = (X' + 1/rho*Th')*P/(P'*P);
    % Update of Th and rho
    Th = Th + rho*(X - P*Q');
    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: Updating rho');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end     
    rho = min([ita*rho, rhoMax]);
    if debug
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('Stage: All Updated');
        disp(['Target Function now ' ': ' num2str(tarFun(X))]);
        disp(['Constraint Deviation now : ' num2str(cstrFun(X, P, Q))]);
        disp(['Largrangian now : ' num2str(largFun(X, P, Q, Th, rho))]);
        
    end 
    % Terminating Condition Check
    if norm(X_old - X, 'fro')/(norm(X_old, 'fro') + 1e-10) < tol ...
            && norm(P_old - P, 'fro')/(norm(P_old, 'fro') + 1e-10) < tol ...
            && norm(Q_old - Q, 'fro')/(norm(Q_old, 'fro') + 1e-10) < tol
        isConverge = true;
    end
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
end

function X = updateX(X, Y, L, X_, R, B, alpha, rho)
%updateX solves the problem: 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} +
%rho/2*||X - X_||_F^2
%
% Where X_ = P*Q' - Th/rho

D = @(X) X - R*X*B;
DY = D(Y);
L_ = eye(size(L)) + 2*alpha*L;
H = DY - R'*DY*B' + rho*X_;

tarFun = @(X) 1/2*norm(DY - D(X), 'fro')^2 + alpha*trace(D(X)'*L*D(X)) + rho/2*norm(X - X_, 'fro');
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
    disp('X unchanged due to non-decreasing within tolerance.');
end
end