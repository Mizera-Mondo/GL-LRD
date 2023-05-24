n = 30;
T = 200;
Y = randn(n, T);
R = eye(n);
k = 30;
alpha = 1; beta = 1;
alpha = alpha/T;
[X, L] = GL_NO_LRD(Y, R, k - 1, alpha, beta);
disp('================================================================');
disp('The graph learning without low-rank decomposition has converged.');
disp('================================================================');

[Xd, Ld] = GL_LRD(Y, R, k - 1, alpha, beta);