function q116(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];


% S0 d f theta phi
lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 

globTol = 0.1; % 0.1 recommended
sigma = eye(5);
sigma(1,1) = 1000;
sigma(2,2) = 0.001;
sigma(3,3) = 1;
sigma(4,4) = pi;
sigma(5,5) = pi;


fmincon_options = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off','GradObj','on','GradConstr','on');%, 'DerivativeCheck', 'on');


model = 'BallStickSSDGrad';

nr_iterations = 10;
% Now run the fitting
% plots the computation time
tic
[parameter_hat, RESNOM] = fitVoxGlobConGeneric(Avox, qhat, bvals, nr_iterations, startx, lb, ub,sigma,globTol, fmincon_options, model)
toc

end