Ita = linspace(0.01, 1, 30);
A = randn(50, 30);
k = 20;
tic
[Q, S, V] = svd(A);
A_svd = Q(:, 1:k)*S(1:k, 1:k)*(V(:, 1:k))';
t_svd = toc;
err_svd = norm(A - A_svd, 'fro');
Err = zeros(size(Ita));
T = Err;
parfor i = 1:30
    disp(['Starting ' num2str(i) '-th iteration']);
    tic
    [P, Q] = LRD_normed(A, k, ita = Ita(i), solver = 'closedform');
    T(i) = toc;
    Err(i) = norm(A - P*Q', 'fro');
end
close all
semilogy(Ita, T);
hold on
semilogy(Ita, t_svd*ones(1, 30));
grid
title('Time');
legend('Normed PQ', 'SVD Baseline');
xlabel('ita');
ylabel('s');
figure;
semilogy(Ita, Err.^2);
hold on
semilogy(Ita, err_svd.^2*ones(1, 30));
ylim([err_svd.^2-1, max(Err.^2)]);
xlabel('ita');
grid
title('NMSE');
legend('Normed PQ', 'SVD Baseline');