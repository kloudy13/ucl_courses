function [two_sigma, conf95] = q121()

[dwis, qhat, bvals] = q1Preprocessing();

NR_PARAMS = 3;
Avox = dwis(:,52,62,25); %given voxel
%Avox = dwis(:,63,40,18);
%Avox = dwis(:,70,64,14);
%Avox = dwis(:,111,111,14);
%Avox = dwis(:,99,63,25);
%Avox = dwis(:,32,15,25);

sample_data(dwis, qhat, bvals, Avox);

%load('q121.mat');
load('q121.mat');

two_sigma = zeros(NR_PARAMS, 2);
conf95 = zeros(NR_PARAMS, 2);

for p=1:NR_PARAMS
    [h, two_sigma(p,:), conf95(p, :)] = q12calcUncertainty(parameter_sample(:,p), 1);
    filename = sprintf('report/figures/q2/q121-p%d.eps', p);
    %hgexport(h, filename);
end

end

function sample_data(dwis, qhat, bvals, Avox)

[NR_IMAGES, ~, ~, ~] = size(dwis);

nr_iterations = 10;

startx = [1.1e+05 2e-03 0.5 0 0];
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

[parameter_hat, minSSD] = fitVoxGlobUnc(Avox, qhat, bvals, nr_iterations, startx, sigma, fminuncOptions, globTol, model);
predicted = BallStick(parameter_hat, bvals, qhat);

T = 1000;
NR_PARAMS = length(parameter_hat);
sigma = sqrt(sum((Avox - predicted').^2)/(NR_IMAGES-NR_PARAMS));

% Define various options for the non-linear fitting algorithm
options = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter', 'Display', 'off');
% bounds for S0 d f theta phi 
lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 

parameter_sample = zeros(T,5);
for t=1:T
  t
  Ahat = predicted + normrnd(0, sigma, size(predicted));
  % Now run the fitting
  parameter_sample(t,:) = fmincon('BallStickSSD', startx, [],[],[],[],lb, ub, [], options, Ahat', bvals, qhat);
end

save('q121.mat', 'parameter_sample');


end

