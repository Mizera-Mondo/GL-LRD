Ita = 0.03;
signal_length = 10000;
node_num = 20:20:1000;
A = randn(10000, 100);

Err = zeros(size(node_num));
T = Err;
Err_svd = Err;
T_svd = T;
for i = 1:length(node_num)
    disp(['Starting ' num2str(i) '-th iteration']);
    A = randn(node_num(i), signal_length);
    k = ceil(0.5*node_num(i));

    tic
    [Q, S, V] = svd(A);
    A_svd = Q(:, 1:k)*S(1:k, 1:k)*(V(:, 1:k))';
    T_svd(i) = toc;
    Err_svd(i) = norm(A - A_svd, 'fro');

    tic
    [P, Q] = LRD_normed(A, k, ita = Ita, solver = 'closedform');
    T(i) = toc;
    Err(i) = norm(A - P*Q', 'fro');

end
save("results_size.mat");

close all
semilogy(node_num, T, LineWidth = 1);
hold on
semilogy(node_num, T_svd, LineWidth = 1);
ylim([min([T_svd, T], [], 'all')./1.5, max([T_svd, T], [], 'all').*1.5]);
grid
title('Time');
legend('Alternative', 'SVD');
xlabel('Number of Nodes');
ylabel('Time /s');

figure;
semilogy(node_num, Err.^2, LineWidth = 1.5);
hold on
semilogy(node_num(1:2:end), (Err_svd(1:2:end)).^2,'+', LineWidth = 2);
xlabel('Number of Nodes');
ylim([min([Err_svd, Err], [], 'all').^2/1.5, max([Err_svd, Err], [], 'all').^2*1.5]);
yyaxis right
relative_err = 100*abs(Err.^2 - Err_svd.^2)./Err_svd.^2;
bar(node_num(1:2:end), relative_err(1:2:end));
ylim([0, 2*max(relative_err)]);
grid
title('NMSE');
legend('Alternative', 'SVD', 'Relative Diff %');