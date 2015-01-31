function [sumRes, S] = q13ZeppelinStickTortTransSSD(x, Avox, bvals, qhat)

[S0, diff, f, theta, phi, lam1, lam2] = q3ZeppTrans(x);

% Synthesize the signals
fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGradSquared = (sum(qhat .* repmat(fibDir, [length(qhat) 1])')).^2;

Si = exp(-bvals * diff .* fibDotGradSquared); % intra-cellular diffusion
Se = exp(-bvals.*(lam2 + (lam1 - lam2)*fibDotGradSquared)); % extra-cellular diffusion

S = S0*(f*Si + (1-f)*Se);

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

end