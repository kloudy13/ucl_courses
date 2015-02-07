function q116(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

nr_iterations = 10;
% Now run the fitting
% plots the computation time
tic
[parameter_hat, RESNOM, minCounter] = fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx)
toc


end
