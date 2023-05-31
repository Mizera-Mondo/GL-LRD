close all
examNum = 5;
nodeNum = 30;
usedEigNum = 5;
signalLength = 100;
noiseCov = 0;
rPerturbation = 0.01;
demoSignalFeature(nodeNum, examNum, usedEigNum, signalLength, noiseCov, rPerturbation);
examNum = 5;
nodeNum = 30;
usedEigNum = 5;
signalLength = 100;
noiseCov = 1;
rPerturbation = 0.01;
demoSignalFeature(nodeNum, examNum, usedEigNum, signalLength, noiseCov, rPerturbation);