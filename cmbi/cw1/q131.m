function q131()

[signals, bvals, qhat] = q13preprocessing();


% Define a starting point for the non-linear fit
%startx = [1.1e+05 1.5e-03 0.5 0 0];

% fit the linear model first
[parameter_hat1, minSSD1] = q1FitVoxLin(signals, qhat, bvals)

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-15, 'TolFun', 1e-15, 'Display', 'iter');
% S0 d f theta phi

nr_iterations = 1;
% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
tic
[parameter_hat2, minSSD2] = q3fitVoxGlobCon(signals, qhat, bvals, nr_iterations, parameter_hat1)
toc
% plots the computation time

end