function [parameter_hat2, minSSD2] = q131()

[signals, bvals, qhat] = q13preprocessing();


% Define a starting point for the non-linear fit
%startx = [1.1e+05 1.5e-03 0.5 0 0];

% fit the linear model first
%[parameter_hat1, minSSD1] = q1FitVoxLin(signals, qhat, bvals)

% the following starting point was obtained after doing linear fit followed by 
% repeated calls to q3fitVoxGlobCon, minSSD is around 4.536
parameter_hat1 = [0.897793040801267   0.000000001010154   0.460722376201871   4.731145749356769, 0.018281659676733];

h = eyeball(signals, parameter_hat1, bvals, qhat);

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-100, 'TolFun', 1e-100, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 5;
% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
tic
[paramsBallStick, SSDDBallStick] = q3fitVoxGlobCon(signals, qhat, bvals, nr_iterations, parameter_hat1)
toc
% plots the computation time

h = eyeball(signals, paramsBallStick, bvals, qhat);

save('q131BallStick.mat', 'paramsBallStick', 'SSDDBallStick');

end