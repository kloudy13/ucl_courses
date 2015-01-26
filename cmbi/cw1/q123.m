function q123()

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

startx = [1.1e+05 2e-03 0.5 0 0];
nr_iterations = 20;
twoSigmaLap = LaplaceUncert(Avox, qhat, bvals, nr_iterations, startx)';

[two_sigmaParBoot, conf95ParBoot] = q121();
[two_sigmaMCMC, conf95MCMC] = q122();

for p=1:NR_PARAMS
    plot(twoSigmaLap(:,p), [1,1]);
    hold on
    plot(two_sigmaParBoot(:,p), [2,2]);
    hold on
    plot(two_sigmaMCMC(:,p), [3,3]);
    legend('2sigma Laplace','2 sigma Par Bootstrap', '2sigma MCMC');
end

%  diag cov - 2sigma range
%   -0.166775137824368
%   -0.000000000000000
%   -0.000000000016855
%   -0.000000000018668
%   -0.000000000025459

end

function twoSigma = LaplaceUncert(Avox, qhat, bvals, nr_iterations, startx)

[parameter_hat, ~, ~, Hessian] = fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx);

cov = -inv(Hessian);

twoSigma = [parameter_hat - diag(cov(1)); parameter_hat + diag(cov(1))];
twoSigma = twoSigma(:,1:3);

end