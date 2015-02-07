function q113(dwis, qhat, bvals)
[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25); %given voxel
%Avox = dwis(:,63,40,18);
%Avox = dwis(:,50,64,23);
%Avox = dwis(:,111,111,14);
%Avox = dwis(:,99,63,25);
%Avox = dwis(:,32,15,25);

nr_iterations = 100;

startx = [1.5 3e-03 0.5 0 0];
%startx = [1.132e+05 1.534e-03 0.57 0 0];

% tried tol of 1e-08 to se if we still get errors
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
globTol = 0.1; % 0.1 recommended

sigma = eye(5);
sigma(1,1) = 0.5 * sqrt(1.1e+05);
sigma(2,2) = 0.2 * sqrt(1e-03);
sigma(3,3) = 1;
sigma(4,4) = sigAngleScale *pi;
sigma(5,5) = sigAngleScale *pi;
sigma = sigma .^2;
model = 'BallStickSSDq112';

tic
[parameter_hat, minSSD, minCounter] = fitVoxGlobUnc(Avox, qhat, bvals, nr_iterations, startx, sigma, fminuncOptions, globTol, model)
toc

predicted = BallStick(parameter_hat, bvals, qhat);
h = eyeball(Avox, predicted, bvals, qhat);
hgexport(h, 'report/figures/q1/q113.eps');

parameter_hat(1)
parameter_hat(2)
parameter_hat(3)
parameter_hat(4)
parameter_hat(5)
end