function [sumRes, S] = q118BallStickRician(x, Avox, bvals, qhat, sigma)

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

% Compute the Riccian noise model

I0 = log(besseli(0,(Avox' .* S)/sigma^2));
sumRes = sum(-(Avox'.^2 + S.^2)/(2*sigma^2)) + sum(I0);

%sumRes = abs(sumRes); % convert from complex to real
end