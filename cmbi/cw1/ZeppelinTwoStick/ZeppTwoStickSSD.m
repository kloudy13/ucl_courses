function [sumRes, S] = ZeppTwoStickSSD(x, Avox, bvals, qhat)

% Extract the parameters
S0 = x(1);
diff = x(2);
f1 = x(3);
f2 = x(4);
theta1 = x(5);
phi1 = x(6);
theta2 = x(7);
phi2 = x(8);
lam1 = x(9);
lam2 = x(10);

% Synthesize the signals
fibDir1 = [cos(phi1)*sin(theta1) sin(phi1)*sin(theta1) cos(theta1)];
fibDir2 = [cos(phi2)*sin(theta2) sin(phi2)*sin(theta2) cos(theta2)];

fibDotGradSquared1 = (sum(qhat .* repmat(fibDir1, [length(qhat) 1])')).^2;
fibDotGradSquared2 = (sum(qhat .* repmat(fibDir2, [length(qhat) 1])')).^2;

Si1 = exp(-bvals * diff .* fibDotGradSquared1); % intra-cellular diffusion
Si2 = exp(-bvals * diff .* fibDotGradSquared2); % intra-cellular diffusion

f1e = f1/(f1+f2); % influence of first stick, normalised 
f2e = f2/(f1+f2); % influence of second stick, normalised
Se = exp(-bvals.*(lam2 + (lam1 - lam2)*(f1e*fibDotGradSquared1+f2e*fibDotGradSquared2))); % extra-cellular diffusion

S = S0*(f1*Si1 + f2*Si2 + (1-f1-f2)*Se);

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

end