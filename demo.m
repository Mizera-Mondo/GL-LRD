n = 30;
T = 200;
Y = randn(n, T);
R = eye(n);
k = 30;
alpha = 1; beta = 1;
alpha = alpha/T;
[X, L] = GL_LRD_NORM_PQ(Y, R, k, alpha, beta);