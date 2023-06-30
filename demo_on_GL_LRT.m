nodeNum = 30;
usedEigNum = 25;
signalLength = 1000;
noiseCov = 0.1;
rPerturbation = 0.05;
[Y, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPerturbation);
[X, Lest_, Aest] = GL_LRT(Y, R, usedEigNum, alpha = 1, beta = 5, LowRankApprox = true);
close all
imagesc(A); title('Ground Truth');
figure; 
imagesc(Aest); title('Estimated');
figure;
[~, Z, ~] = svd(Y); imagesc(Z(:, 1:nodeNum));  title('Singular Value of Y');
figure; 
[~, Z, ~] = svd(X); imagesc(Z(:, 1:nodeNum)); title('Singular Value of X');