function [parameter_hat] = mapDiffTenToBallStick(D, logS0)
% D is a 3x3 matrix with diffusivity in X, Y, Z directions

[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
[maxEig, maxEigPos] = max(EigVals);

principalDir = EigVect(:,maxEigPos);

%theta = acos(principalDir(3));
%phi = acos(principalDir(1)/sin(theta));

[phi, theta]= cart2sph(principalDir(1),principalDir(2),principalDir(3))
S0 = exp(logS0);

d = mean(D(:));

EigValsScaled = EigVals/ sum(EigVals);

f = 0.5 * (abs(EigValsScaled(1) - EigValsScaled(2))...
         + abs(EigValsScaled(1) - EigValsScaled(3))...
         + abs(EigValsScaled(2) - EigValsScaled(3)));

parameter_hat = [S0 d f theta phi];