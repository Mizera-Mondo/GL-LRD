n = 30;
T = 200;
Y = randn(n, T);
R = eye(n);
k = 10;
alpha = 1; beta = 1;
alpha = alpha/T;
[X, L] = LR_DGI(Y, R, k, alpha, beta);