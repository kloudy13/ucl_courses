function q116(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi
lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 

nr_iterations = 10;
% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
tic
[parameter_hat, RESNOM] = fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx)
toc
% plots the computation time

end