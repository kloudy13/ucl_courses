function q112(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();
Avox = dwis(:,52,62,25);

% apply the inverse tranformations: sqrt and tangent
startx = [sqrt(7.5e+05) sqrt(3e-03) q1TransInv(2.5e-01) 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSDq112', startx, h, Avox, bvals, qhat);
RESNOM

% apply the transformations
[S0 d f theta phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
parameter_hat = [ S0^2 d^2 q1Trans(f) theta phi]

h = eyeball(Avox, parameter_hat, bvals, qhat);


end