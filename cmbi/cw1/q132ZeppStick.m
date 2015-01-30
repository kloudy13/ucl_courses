function [paramsZeppStick, SSDZeppStick] = q132ZeppStick()

[signals, bvals, qhat] = q13preprocessing();

params_orig = [0.897793040801267, 0.000000001010154, 0.460722376201871, 4.731145749356769, 0.018281659676733, 6.73436e-10, 1.346872e-09];
[ ~, predicted ]= q13ZeppelinStickSSD(params_orig, signals, bvals, qhat);
%h = eyeball(signals, predicted, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-100, 'TolFun', 1e-100, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 150;
lb = [0  , 0  , 0, -inf, -inf, 0, 0];
ub = [inf, inf, 1,  inf,  inf, inf, inf]; 
fminconOptions = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');

sigma = eye(5);
sigma(1,1) = 0.15;
sigma(2,2) = 1e-09;
sigma(3,3) = 0.15;
sigma(4,4) = pi;
sigma(5,5) = pi;
sigma(6,6) = 1e-09;
sigma(7,7) = 1e-09;
model = 'q13ZeppelinStickSSD';

% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
tic
[paramsZeppStick, SSDZeppStick] = q3fitVoxGlobCon(signals, qhat, bvals, nr_iterations, params_orig, lb, ub, sigma, fminconOptions, model)
toc
% plots the computation time

[ ~, predicted ]= q13ZeppelinStickSSD(paramsZeppStick, signals, bvals, qhat);
h = eyeball(signals, predicted, bvals, qhat);

save('q132DiffTensor.mat', 'paramsZeppStick', 'SSDZeppStick');
params_orig - paramsZeppStick
end