nodeNum = 10;
usedEigNum = 10;
signalLength = 1000;
noiseCov = 0;
rPerturbation = 0.05;
[Y, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPerturbation);
[X, Lest, Aest] = GL_LRT(Y, R, usedEigNum, 1, 2, LowRankApprox = false);
close all
imagesc(A); figure; imagesc(Aest);