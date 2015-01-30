function [logS0, D, SSDDiffTensor] = q132DiffTensor()

[signals, bvals, qhat] = q13preprocessing();

[logS0, D, SSDDiffTensor] = q1FitVoxLin(signals, qhat, bvals);

save('q132DiffTensor.mat', 'logS0', 'D', 'SSDDiffTensor');

end