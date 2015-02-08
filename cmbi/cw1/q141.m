function q141()

[signals, bvals, qhat] = q13preprocessing();

load('q131BallStick.mat');
sigma = 0.04;
F = Fisher(paramsBallStick, signals, bvals, qhat, sigma);

end

function F = Fisher(paramsBallStick, signals, bvals, qhat, sigma)

dSdS0 = BallStickDerivSigS0(paramsBallStick, bvals, qhat);
dSdD = BallStickDerivSigD(paramsBallStick, bvals, qhat);
dSdF = BallStickDerivSigF(paramsBallStick, bvals, qhat);
dSdTheta = BallStickDerivSigTheta(paramsBallStick, bvals, qhat);
dSdPhi = BallStickDerivSigPhi(paramsBallStick, bvals, qhat);

dS = [dSdS0; dSdD; dSdF; dSdTheta; dSdPhi]';

NR_IMAGES = length(signals);
NR_PARAMS = length(paramsBallStick);
F = zeros(NR_PARAMS, NR_PARAMS);
for k=1:NR_IMAGES
  F = F + dS(k,:)' * dS(k,:);
end

F = F / sigma^2;
F = F ./ (paramsBallStick' * paramsBallStick);

end