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

T = 10;
sigma = sum((Avox - predicted).^2)/(NR_IMAGES-length(predicted));

parameter_sample = zeros(T,5);
for t=1:T
  Ahat = predicted + normrnd(0, sigma, size(predicted));
  parameter_sample(t,:) = %use fminunc directly (Ahat, qhat, bvals, nr_iterations, startx);
end

std(parameter_sample)
end