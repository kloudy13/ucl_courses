function q115(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

%voxSig = dwis(:,52,62,25);
voxSig = dwis(:,63,40,18);
%voxSig = dwis(:,70,64,14);
%voxSig = dwis(:,111,111,14);
%voxSig = dwis(:,99,63,25);
%voxSig = dwis(:,32,15,25);

[logS0, D, minSSDLin] = q1FitVoxLin(voxSig, qhat, bvals); 
[param_start] = mapDiffTenToBallStick(D, logS0);
minSSDstart = BallStickSSD(param_start, voxSig, bvals, qhat);

nr_iterations = 100;
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
globTol = 0.1; % 0.1 recommended

sigma = eye(5);
sigma(1,1) = sigmaScale*sqrt(7.5e+05);
sigma(2,2) = sigmaScale*sqrt(3e-03);
sigma(3,3) = sigmaScale*5;
sigma(4,4) = sigAngleScale *pi;
sigma(5,5) = sigAngleScale *pi;
model = 'BallStickSSDq112';

[parameter_hat2, minSSD2, minCounter] = fitVoxGlobUnc(voxSig, qhat, bvals,nr_iterations, param_start, sigma, fminuncOptions, globTol, model);

minCounter
param_start - parameter_hat2
minSSDstart - minSSD2


end