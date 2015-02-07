function q115(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

voxSig = dwis(:,52,62,25);
%voxSig = dwis(:,63,40,18);
%voxSig = dwis(:,70,64,14);
%voxSig = dwis(:,111,111,14);
%voxSig = dwis(:,99,63,25);
%voxSig = dwis(:,32,15,25);

[logS0, D, minSSDLin] = q1FitVoxLin(voxSig, qhat, bvals); 
% 88, 89,90 times global min was found
%[param_start] = mapDiffTenToBallStick2(D, logS0); % mean of diag(D)

% 94,92,92 times global min was found
[param_start] = mapDiffTenToBallStick(D, logS0); %mean of D(:)

% 86,90,89 times global min was found
%[param_start] = mapDiffTenToBallStick3(D, logS0); %max(D(:)), f uses (lam1-lam2)^2, etc ..

minSSDstart = BallStickSSD(param_start, voxSig, bvals, qhat);
%param_start = [1.5 3e-03 0.5 0 0];

nr_iterations = 100;
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
globTol = 0.1; % 0.1 recommended

sigma = eye(5);
sigma(1,1) = 0.05 * sqrt(200);
sigma(2,2) = 0.02 * sqrt(1e-05);
sigma(3,3) = 0.01;
sigma(4,4) = sigAngleScale *pi;
sigma(5,5) = sigAngleScale *pi;
sigma = sigma .^2;
model = 'BallStickSSDq112';

[parameter_hat2, minSSD2, minCounter] = fitVoxGlobUnc(voxSig, qhat, bvals,nr_iterations, param_start, sigma, fminuncOptions, globTol, model)

minCounter



end