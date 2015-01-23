function q121(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

SLICE_NR=25;

[NR_IMAGES,W,H,~] = size(dwis);

Avox = dwis(:,52,62,25); %given voxel
%Avox = dwis(:,63,40,18);
%Avox = dwis(:,70,64,14);
%Avox = dwis(:,111,111,14);
%Avox = dwis(:,99,63,25);
%Avox = dwis(:,32,15,25);

nr_iterations = 10;

startx = [1.1e+05 2e-03 0.5 0 0];
[parameter_hat, minSSD] = fitVoxGlob1(Avox, qhat, bvals, nr_iterations, startx)
predicted = BallStick(parameter_hat, bvals, qhat);

T = 100;
NR_PARAMS = length(parameter_hat);
sigma = sqrt(sum((Avox - predicted').^2)/(NR_IMAGES-NR_PARAMS));
options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 


parameter_sample = zeros(T,5);
for t=1:T
  Ahat = predicted + normrnd(0, sigma, size(predicted));
  need to apply tranformations before calling fminunc
  parameter_sample(t,:) = fminunc('BallStickSSDq112', startx, options, Ahat', bvals, qhat);
  h = eyeball(Ahat, parameter_sample(t,:), bvals, qhat);
end

for p=1:NR_PARAMS
    calcUncertainty(parameter_sample(:,p));    
end

end

function calcUncertainty(parameter_sample)
% input is a 1D parameter vector, respresenting samples for one of the 5
% parameters

mx = mean(parameter_sample);
sx = std(parameter_sample);
N = length(parameter_sample);

% Approximative
upper_limit=mx+1.66*sx;
lower_limit=mx-1.66*sx;
%plot(1:N,parameter_sample,'o',1:N,upper_limit,1:N,lower_limit)

hist(parameter_sample, 20)
end