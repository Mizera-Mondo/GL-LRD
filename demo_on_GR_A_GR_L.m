% 2 different approaches of graph refinement.
% In GL_LRSS, the target variable is L and in GL_LRT, A.
% This demo will turn off low-rank components estimation to compare between
% these 2 approaches.
nodeNum = 20;
usedEigNum = 20;
signalLength = 2000;
noiseCov = 0.01;
rPertubation = 0.01;
[Y, A, R] = genRandomSignal(nodeNum, usedEigNum, signalLength, noiseCov, rPerturbation);
[~, Lest_L] = GL_LRSS(Y, R = R, alpha = 0.3, beta = 5, lowRankEst = false);
[~, Lest_A, Aest_A] = GL_LRT(Y, R, usedEigNum, alpha = 1.5, beta = 3, LowRankApprox = false);
close all
imagesc(diag(sum(A)) - A); title('Ground Truth'); 
figure; imagesc(Lest_L); title('Estimated with GR_L');
figure; imagesc(Lest_A); title('Estimated with GR_A');