function [parameter_hat] = mapDiffTenToBallStick3(D, logS0)
% D is a 3x3 matrix with diffusivity in X, Y, Z directions
% d is calculated as the mean of the diagonal elements of D
% f is calculated using a linear formula, being 0 if all the eigenvalues are the
% same and 1 if 2 eigenvalues are zero. In between it interpolates linearly.


[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
[maxEig, maxEigPos] = max(EigVals);

principalDir = EigVect(:,maxEigPos);

%theta = acos(principalDir(3));
%phi = acos(principalDir(1)/sin(theta));

[phi, theta]= cart2sph(principalDir(1),principalDir(2),principalDir(3))
S0 = exp(logS0);

d = mean(max(D(:)));

EigValsScaled = EigVals/ sum(EigVals);

f = 0.5 * ((EigValsScaled(1) - EigValsScaled(2))^2 ...
         + (EigValsScaled(1) - EigValsScaled(3))^2 ...
         + (EigValsScaled(2) - EigValsScaled(3))^2);

parameter_hat = [S0 d f theta phi];