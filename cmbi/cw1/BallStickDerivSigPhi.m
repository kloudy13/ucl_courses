function [dSdPhi] = BallStickDerivSigPhi(x, bvals, qhat)

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGrad = sum(qhat .* repmat(fibDir, [length(qhat) 1])');

Si = exp(-bvals * diff .* (fibDotGrad .^2));

A = S0*f*Si .*(-(2*bvals*diff) .* fibDotGrad);

derivQN = repmat([-sin(phi)*sin(theta) cos(phi)*sin(theta) 0], [length(qhat) 1]);

dSdPhi = A .* sum(derivQN' .* qhat);

end