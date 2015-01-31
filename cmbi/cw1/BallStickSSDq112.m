function [sumRes, S] = BallStickSSDq112(x, Avox, bvals, qhat)
% same as BallStickSSD but with transformations of input parameters in 
% order to satisfy boundary constraints 

% Extract the parameters
% make S0 and diff positive
S0 = x(1)^2;
diff = x(2)^2;

% f must be between 0,1 .. so apply the inverse tangent function
%f = atan(x(3))/pi+0.5;
f = q1Trans(x(3));
%f = x(3);

% leave theta and phi unconstrained
theta = x(4);
phi = x(5);

% Synthesize the signals
fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGrad = sum(qhat .* repmat(fibDir, [length(qhat) 1])');

S = S0*(f*exp(-bvals * diff .* (fibDotGrad .^2)) + (1-f)*exp(-bvals*diff));

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

end