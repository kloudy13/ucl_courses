function [paramsNoddiSimple, SSDNoddiSimple] = NoddiSimpleUnc()

addpath ../

[signals, bvals, qhat] = q13preprocessing();

[signals, bvals, qhat] = q13preprocessing();
[logS0, D, SSDDiffTensor] = q1FitVoxLin(signals, qhat, bvals);
[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
sortedEigs = sort(EigVals);

% S0 = 1;
% d = 1e-09;
% f = 0.34;
% theta = 4.68;
% phi = 3.166;
% k = 0.5;

S0 = 0.923198831615024;
d = 0.000000001135595 ;
f = 0.499421794318001;
theta = 4.704686705926314;
phi = 0.038495126169283;
k = 15;

params_orig = [S0, d, f, theta, phi, k];

%load('NoddiSimple.mat');
%params_orig = paramsNoddiSimple;

[ orig_ssd, predicted ]= NoddiSimpleSSD(params_orig, signals, bvals, qhat);
h = eyeball(signals, predicted, bvals, qhat);
orig_ssd
%optimal_params = [0.9816, 6.339e-10, 0.3464, 1.5468, -3.1183, 2.8719e-09, 6.868e-10]
%[ ~, predicted_opt ]= NoddiSimpleSSD(optimal_params, signals, bvals, qhat);
%h = eyeZepp(signals, predicted_opt, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 2000, 'Algorithm', 'interior-point',...
    'TolX', 1e-20, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 0;
globTol = 0.1; % 0.1 recommended

sigmaScale = 0.1;

sigma = eye(5);
sigma(1,1) = 0.15; %S0
sigma(2,2) = 1e-09; %d
sigma(3,3) = 0.1; %f
sigma(4,4) = 2*pi; %theta
sigma(5,5) = 2*pi; %phi
sigma(6,6) = 1; %k
%sigma(9,9) = 1e-09; %lam1
%sigma(10,10) = 1e-09; %lam2

model = 'NoddiSimpleTransSSD';
trans = str2func('NoddiSimpleTrans');
transInv = str2func('NoddiSimpleTransInv');
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter'); 

% Now run the fitting
tic
[paramsNoddiSimple, SSDNoddiSimple] = q3fitVoxGlobUnc(signals, qhat, bvals, nr_iterations, params_orig, sigma, fminuncOptions, globTol, model, trans, transInv)
toc
% plots the computation time

[ ~, predicted ]= NoddiSimpleSSD(paramsNoddiSimple, signals, bvals, qhat);
h = eyeball(signals, predicted, bvals, qhat);

save('NoddiSimple.mat', 'paramsNoddiSimple', 'SSDNoddiSimple');
% best SSD is 1.1534
end