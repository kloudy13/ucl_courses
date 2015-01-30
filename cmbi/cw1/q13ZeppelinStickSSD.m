function [sumRes, S] = q13ZeppelinStickSSD(x, Avox, bvals, qhat)

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
lam1 = x(6); % lam1 > lam2
lam2 = x(7);

% Synthesize the signals
fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGradSquared = (sum(qhat .* repmat(fibDir, [length(qhat) 1])')).^2;

Si = exp(-bvals * diff .* fibDotGradSquared); % intra-cellular diffusion
Se = exp(-bvals.*(lam2 + (lam1 - lam2)*fibDotGradSquared)); % extra-cellular diffusion

S = S0*(f*Si + (1-f)*Se);

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

end