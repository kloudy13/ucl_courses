function pAgX = BallStickProb(x, Avox, bvals, qhat, sigma)
% computes p(x|A), where x is the predicted model and A are the measurements

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
pAgX = prod(normpdf(S, Avox, sigma * ones(size(S))));

end