function [paramsZeppTwoStick, SSDZeppTwoStick] = ZeppTwoStickUnc()

addpath ../

[signals, bvals, qhat] = q13preprocessing();

[signals, bvals, qhat] = q13preprocessing();
[logS0, D, SSDDiffTensor] = q1FitVoxLin(signals, qhat, bvals);
[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
sortedEigs = sort(EigVals);

S0 = 1;
d = 1e-09;
f1 = 0.3;
f2 = 0.2;
theta1 = 4;
phi1 = 0.02;
theta2 = 2;
phi2 = -1;
lam1 = sortedEigs(3);
lam2 = sortedEigs(2);


%params_orig = [S0, d, f1, f2, theta1, phi1, theta2, phi2, lam1, lam2];

load('ZeppTwoStick.mat');
params_orig = paramsZeppTwoStick;

%[ ~, predicted ]= ZeppTwoStickSSD(params_orig, signals, bvals, qhat);
%h = eyeZepp(signals, predicted, bvals, qhat);

%optimal_params = [0.9816, 6.339e-10, 0.3464, 1.5468, -3.1183, 2.8719e-09, 6.868e-10]
%[ ~, predicted_opt ]= ZeppTwoStickSSD(optimal_params, signals, bvals, qhat);
%h = eyeZepp(signals, predicted_opt, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 100;
globTol = 0.1; % 0.1 recommended

% sigma = eye(5);
% sigma(1,1) = 0.15; %S0
% sigma(2,2) = 1e-09; %d
% sigma(3,3) = 0.1; %f1
% sigma(4,4) = 0.1; %f2
% sigma(5,5) = 2*pi; %theta1
% sigma(6,6) = 2*pi; %phi1
% sigma(7,7) = 2*pi; %theta2
% sigma(8,8) = 2*pi; %phi2
% sigma(9,9) = 1e-09; %lam1
% sigma(10,10) = 1e-09; %lam2

sigmaScale = 0.1;

sigma = eye(5);
sigma(1,1) = sigmaScale * 0.15; %S0
sigma(2,2) = sigmaScale * 1e-09; %d
sigma(3,3) = sigmaScale * 0.1; %f1
sigma(4,4) = 0.1; %f2
sigma(5,5) = sigmaScale * 2*pi; %theta1
sigma(6,6) = sigmaScale * 2*pi; %phi1
sigma(7,7) = 2*pi; %theta2
sigma(8,8) = 2*pi; %phi2
sigma(9,9) = sigmaScale * 1e-09; %lam1
sigma(10,10) = sigmaScale * 1e-09; %lam2

model = 'ZeppTwoStickTransSSD';
trans = str2func('ZeppTwoStickTrans');
transInv = str2func('ZeppTwoStickTransInv');
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

% Now run the fitting
tic
[paramsZeppTwoStick, SSDZeppTwoStick] = q3fitVoxGlobUnc(signals, qhat, bvals, nr_iterations, params_orig, sigma, fminuncOptions, globTol, model, trans, transInv)
toc
% plots the computation time

[ ~, predicted ]= ZeppTwoStickSSD(paramsZeppTwoStick, signals, bvals, qhat);
h = eyeball(signals, predicted, bvals, qhat);

save('ZeppTwoStick.mat', 'paramsZeppTwoStick', 'SSDZeppTwoStick');
% best SSD is 1.1534
end