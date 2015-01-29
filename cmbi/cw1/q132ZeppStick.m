function [paramsZeppStick, SSDZeppStick] = q132ZeppStick()

[signals, bvals, qhat] = q13preprocessing();

parameter_hat1 = [0.897793040801267, 0.000000001010154, 0.460722376201871, 4.731145749356769, 0.018281659676733, 6.73436e-10, 1.346872e-09];

h = eyeball(signals, parameter_hat1, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-100, 'TolFun', 1e-100, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 5;
lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 
fminconOptions = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');

sigma = eye(5);
sigma(1,1) = 0.15;
sigma(2,2) = 1e-09;
sigma(3,3) = 0.15;
sigma(4,4) = pi;
sigma(5,5) = pi;
model = 'BallStickSSD';

% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
tic
[paramsZeppStick, SSDZeppStick] = q3fitVoxGlobCon(signals, qhat, bvals, nr_iterations, parameter_hat1, lb, ub, sigma, fminconOptions, model)
toc
% plots the computation time

h = eyeball(signals, paramsZeppStick, bvals, qhat);

save('q132DiffTensor.mat', 'paramsZeppStick', 'SSDZeppStick');

end