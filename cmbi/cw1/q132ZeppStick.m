function [paramsZeppStick, SSDZeppStick] = q132ZeppStick()

[signals, bvals, qhat] = q13preprocessing();
continuuee here

[paramsZeppStick, SSDZeppStick] = q1FitVoxLin(signals, qhat, bvals);

save('q132DiffTensor.mat', 'paramsDiffTensor', 'SSDDiffTensor');

end