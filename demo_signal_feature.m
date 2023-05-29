examNum = 5;
nodeNum = 20;
usedEigNum = 5;
signalLength = 100;
sAcc = zeros(nodeNum, 1);
close all
figure
hold on
% yyaxis left
for i = 1:examNum

    A = rand_ugraph(nodeNum, nodeNum*3, 0.1, 0.1);
    A = A./sum(A, "all").*nodeNum;
    L = diag(sum(A)) - A;

    [vec, val] = eig(L);
    [vec, val] = sortEigen(vec, val, 'descend');
    vecT = vec(:, 1:usedEigNum);
    valT = val(1:usedEigNum, 1:usedEigNum);
    covZ = diag(ones(usedEigNum, 1)./(diag(valT)));

    z = randn(usedEigNum, signalLength);
    z = sqrt(covZ)*z;
    v = vecT*z;

    x = zeros(nodeNum, signalLength);
    x(:, 1) = v(:, 1);
    for t = 2:signalLength
        x(:, t) = x(:, t - 1) + v(:, t);
    end

    [~, S, ~] = svd(x);
    s = diag(S);
    sAcc = sAcc + s;
    s(abs(s) < 1e-6) = 0;
    plot(s, '.', 'MarkerSize', 15)
    
end

set(gca, 'xscal', 'log')
% yyaxis right
plot(sAcc./examNum, 'LineWidth', 2)
A = cell(examNum + 1, 1);
for i = 1:examNum
    A{i} = ['Simu. ' num2str(i)];
end
A{i + 1} = 'Average';
legend(A)
title('Singular Spectrum of Signal X')
xlabel('Singular Value Count')
ylabel('Singular Value')
grid on
grid minor
hold off