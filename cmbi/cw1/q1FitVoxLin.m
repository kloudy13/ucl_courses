
function [parameter_hat, minSSD] = q1FitVoxLin(voxSig, qhat, bvals)
% linear model, but doesn't ensure that D is positive definite. Try cholesky:
% D=L*L^T
N = length(voxSig);

[qX, qY, qZ] = deal(qhat(1,:),qhat(2,:),qhat(3,:));

Y = [ones(1,N); -bvals .* qX.^2; -2*bvals.*qX.*qY; -2*bvals.*qX.*qZ; -bvals.*qY.^2; -2*bvals.*qY.*qZ; -bvals.*qZ.^2]';

X = Y\log(voxSig);

[logS0, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz] = deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7));
D = [Dxx Dxy Dxz; Dxy Dyy Dyz; Dxz Dyz Dzz ];

[EigVect,EigVals] = eig(D); 

EigVals = diag(EigVals);
[maxEig, maxEigPos] = max(EigVals);

principalDir = EigVect(:,maxEigPos);

theta = acos(principalDir(3));
phi = acos(principalDir(1)/sin(theta));

S0 = exp(logS0);

d = sum(D(:));

EigValsScaled = EigVals/ sum(EigVals);

f = 0.5 * (abs(EigValsScaled(1) - EigValsScaled(2)) + abs(EigValsScaled(1) - EigValsScaled(3)) + abs(EigValsScaled(2) - EigValsScaled(3)));

parameter_hat = [S0 d f theta phi];

minSSD = BallStickSSD(parameter_hat, voxSig, bvals, qhat);
end
