n = 30;
T = 200;
Y = randn(n, T);
R = eye(n);
k = 30;
[P, Q] = LRD(Y, k);
[Pd, Qd] = LRD(Y, k - 1);