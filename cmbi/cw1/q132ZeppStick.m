function [paramsZeppStick, SSDZeppStick] = q132ZeppStick()

[signals, bvals, qhat] = q13preprocessing();
[logS0, D, SSDDiffTensor] = q1FitVoxLin(signals, qhat, bvals);
[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
sortedEigs = sort(EigVals);
lam1 = sortedEigs(3);
lam2 = sortedEigs(2);

S0 = 1;
d = 1e-09;
f = 0.46;
theta = 4;
phi = 0.02;

params_orig = [S0, d, f, theta, phi, lam1, lam2];
[ ~, predicted ]= q13ZeppelinStickSSD(params_orig, signals, bvals, qhat);
%h = eyeball(signals, predicted, bvals, qhat);
optimal_params = [0.9816, 6.339e-10, 0.3464, 1.5468, -3.1183, 2.8719e-09, 6.868e-10]
[ ~, predicted_opt ]= q13ZeppelinStickSSD(optimal_params, signals, bvals, qhat);
h = eyeball(signals, predicted_opt, bvals, qhat);

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
sigma(3,3) = 0.1;
sigma(4,4) = 2*pi;
sigma(5,5) = 2*pi;
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

save('q132ZeppStick.mat', 'paramsZeppStick', 'SSDZeppStick');
params_orig - paramsZeppStick
% result should be 1.17
end