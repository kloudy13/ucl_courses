function q123()

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

NR_PARAMS = 3;
startx = [1.1e+05 2e-03 0.5 0 0];
nr_iterations = 20;
%twoSigmaLap = LaplaceUncert(Avox, qhat, bvals, nr_iterations, startx)';

[two_sigmaParBoot, conf95ParBoot] = q121Conf();
[two_sigmaMCMC, conf95MCMC] = q122Conf();

for p=2:2 %1:NR_PARAMS
    h = figure
    %plot(twoSigmaLap(:,p), [1,1]);
    %hold on
    plot(two_sigmaParBoot(p,:), [2,2]);
    hold on
    plot(two_sigmaMCMC(p,:), [3,3]);
    hold on
    plot(conf95ParBoot(p,:), [4,4]);
    hold on
    plot(conf95MCMC(p,:), [5,5]);
    ylim ([1 6]);
    %legend('2sigma Laplace','2 sigma Par Bootstrap', '2sigma MCMC');
    legend('2 sigma Par Bootstrap', '2sigma MCMC', 'conf95ParBoot', 'conf95MCMC');
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


function [two_sigma, conf95] = q121Conf()

load('q121.mat');

NR_PARAMS = 3;
two_sigma = zeros(NR_PARAMS, 2);
conf95 = zeros(NR_PARAMS, 2);
for p=1:NR_PARAMS
    [h, two_sigma(p,:), conf95(p, :)] = q12calcUncertainty(parameter_sample(:,p), 0);  
end

end

function [two_sigma, conf95] = q122Conf()

load('q122.mat');

NR_PARAMS = 3;
two_sigma = zeros(NR_PARAMS, 2);
conf95 = zeros(NR_PARAMS, 2);
for p=1:NR_PARAMS
    [h, two_sigma(p,:), conf95(p, :)] = q12calcUncertainty(samples(p,:), 0);  
end

end
