function q132()
% fits the following models: BallStick, Diffusion Tensor, Zeppelin-Stick
% and Zeppelin-Stick with tortuosity

[signals, bvals, qhat] = q13preprocessing();

% BallStick model
%[paramsBallStick, SSDDBallStick] = q131();
load('q131BallStick.mat');
predicted = BallStick(paramsBallStick, bvals, qhat);
h = q3eyeball(signals, predicted, bvals, qhat);
hgexport(h, 'report/figures/q3/q131.eps');

% Diffusion Tensor model
[logS0, D, SSDDiffTensor, predicted] = q1FitVoxLin(signals, qhat, bvals);
%load('q132DiffTensor.mat');
h = q3eyeball(signals, predicted, bvals, qhat);
hgexport(h, 'report/figures/q3/q132-DiffTensor.eps');

% Zeppelin-Stick
%[paramsZeppStick, SSDZeppStick] = q132ZeppStickUnc();
load('q132ZeppStickUnc.mat');
[~, predicted] = q13ZeppelinStickSSD(paramsZeppStick, signals, bvals, qhat);
h = q3eyeball(signals, predicted, bvals, qhat);
hgexport(h, 'report/figures/q3/q132-ZeppStick.eps');

% Zeppelin-Stick with tortuosity
%[paramsZeppStickTort, SSDZeppStickTort] = q132ZeppStickTort();
load('q132ZeppStickTortUnc.mat');
[~, predicted] = q13ZeppelinStickTortSSD(paramsZeppStick, signals, bvals, qhat);
h = q3eyeball(signals, predicted, bvals, qhat);
hgexport(h, 'report/figures/q3/q132-ZeppStickTort.eps');



end