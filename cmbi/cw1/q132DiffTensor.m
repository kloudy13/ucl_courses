function [paramsDiffTensor, SSDDiffTensor] = q132DiffTensor()

[signals, bvals, qhat] = q13preprocessing();

[paramsDiffTensor, SSDDiffTensor] = q1FitVoxLin(signals, qhat, bvals);

save('q132DiffTensor.mat', 'paramsDiffTensor', 'SSDDiffTensor');

end