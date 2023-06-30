function [X, L, A] = GL_LRT(Y, R, k, options)
%LR-DGI Solve 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} +
%beta/2||L||_F^2

arguments
    Y, R, k double
    options.alpha = 0.1;
    options.beta = 0.1;
    options.graphRefineMethod = 'quadprog';
    options.LowRankApprox = true;
    options.debug = false;
end
alpha = options.alpha;
beta = options.beta;
debug = options.debug;

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

tol = 1e-2;
iter = 1;
maxIter = 1000;
isConverge = false;
isMaxIter = false;

if debug
    disp('=========================================================');
    disp('Starting Graph Learning. Low-rank decomposition involved.');
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
    M = M./T;
    tic
    A = solveSubA(M, alpha, beta);
    toc
    L = diag(sum(A)) - A;

    % Optimizing X
    % X = solveSubX(Y, L, R, B, alpha, k);
    if options.LowRankApprox
        tic
        X = updateX_SVD(X, D(Y), L, R, B, alpha, k, solver = 'GPM');
        toc
    end
    
    isConverge = norm(L_old - L, 'fro')/norm(L_old, 'fro') < tol ...
            && norm(X_old - X, 'fro')/norm(X_old, 'fro') < tol;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
end

function A = solveSubA(M, alpha, beta, options)
%solveSubA solve the sub-problem ||A||_F^2 + alpha/beta*Tr{AM}
arguments
    M, alpha, beta double
    options.method = 'quadprog'
end

if strcmp(options.method, 'quadprog')
    % DON'T USE PERSISTENT VARIABLES TO AVOID DUPLICATED LABOUR OF
    % CONSTRUCTING TARFUN AND CONSTRS. THIS WILL CAUSE BUG SINCE CALLING
    % THIS FUNCTION FOR ANOTHER TIME WILL NOT CLEAR THESE VARIABLES!!!

    [n, ~] = size(M);
    Aeq = [];
    beq = [];
    Aie = [];
    bie = [];

    % Construct target function for vectorized A
    H = eye(n^2);
    f = alpha/beta*mat2vec(M);

    % Construct constraints for vectorized A

    % Equality part: 1'A1 = n, A = A', Aii = 0
    % 1'A1 = n
    Aeq = ones(1, n^2);
    beq = n;
    % A = A'
    for i = 1:n
        for j = 1:i - 1
            aeq = zeros(n);
            aeq(i, j) = 1;
            aeq(j, i) = -1;
            Aeq = [Aeq; (mat2vec(aeq))'];
            beq = [beq; 0];
        end
    end
    % Aii = 0
    for i = 1:n
        aeq = zeros(n);
        aeq(i, i) = 1;
        Aeq = [Aeq; (mat2vec(aeq))'];
        beq = [beq; 0];
    end

    % Inequality part: Aij >= 0, i ~= j
    Aie = -1*eye(n^2);
    bie = zeros(n^2, 1);

    A = quadprog(H, f, Aie, bie, Aeq, beq, [], [], [], optimoptions('quadprog', 'Display','off'));
    A = vec2mat(A, n);

elseif strcmp(options.method, 'CVX')
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

else
    error('%s is not a vaild solver!', options.method);
end

end

