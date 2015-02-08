function [sumRes, gradient] = BallStickSSDGrad(x, Avox, bvals, qhat)

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGrad = sum(qhat .* repmat(fibDir, [length(qhat) 1])');

S = S0*(f*exp(-bvals * diff .* (fibDotGrad .^2)) + (1-f)*exp(-bvals*diff));

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

% Compute the gradient
m2AmS = -2*(Avox' - S);
gradient = zeros(5,1);
gradient(1) = sum(m2AmS .* BallStickDerivSigS0(x, bvals, qhat));
gradient(2) = sum(m2AmS .* BallStickDerivSigD(x, bvals, qhat));
gradient(3) = sum(m2AmS .* BallStickDerivSigF(x, bvals, qhat));
gradient(4) = sum(m2AmS .* BallStickDerivSigTheta(x, bvals, qhat));
gradient(5) = sum(m2AmS .* BallStickDerivSigPhi(x, bvals, qhat));

end