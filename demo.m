n = 30;
T = 200;
Y = randn(n, T);
R = eye(n);
k = 29;
alpha = 1; beta = 1;
alpha = alpha/T;
[X, L] = GL_LRD(Y, R, k, alpha, beta);