function [paramsZeppStick, SSDZeppStick] = q132ZeppStickUnc()

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
%trans = str2func('q3ZeppTrans');
%params_orig2 = trans(params_orig);

[ ~, predicted ]= q13ZeppelinStickSSD(params_orig, signals, bvals, qhat);
%h = eyeball(signals, predicted, bvals, qhat);

%optimal_params = [0.9816, 6.339e-10, 0.3464, 1.5468, -3.1183, 2.8719e-09, 6.868e-10]
%[ ~, predicted_opt ]= q13ZeppelinStickSSD(optimal_params, signals, bvals, qhat);
%h = eyeball(signals, predicted_opt, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 100;
sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
globTol = 0.1; % 0.1 recommended

sigma = eye(5);
sigma(1,1) = 0.15;
sigma(2,2) = 1e-09;
sigma(3,3) = 0.1;
sigma(4,4) = 2*pi;
sigma(5,5) = 2*pi;
sigma(6,6) = 1e-09;
sigma(7,7) = 1e-09;
%sigma = sigma .^2;
model = 'q13ZeppelinStickTransSSD';
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

% Now run the fitting
tic
[paramsZeppStick, SSDZeppStick] = q3fitVoxGlobUnc(signals, qhat, bvals, nr_iterations, params_orig, sigma, fminuncOptions, globTol, model)
toc
% plots the computation time

[ ~, predicted ]= q13ZeppelinStickSSD(paramsZeppStick, signals, bvals, qhat);
h = eyeball(signals, predicted, bvals, qhat);

save('q132ZeppStickUnc.mat', 'paramsZeppStick', 'SSDZeppStick');
params_orig - paramsZeppStick
% ryan got 1.17
end