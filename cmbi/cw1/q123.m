function q123()

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

NR_PARAMS = 3;
startx = [1.1e+05 1.5e-03 0.5 0 0];
nr_iterations = 20;
sigmaNoise = 6000;
twoSigmaLap = LaplaceUncert(Avox, qhat, bvals, nr_iterations, startx, sigmaNoise);

[two_sigmaParBoot, conf95ParBoot] = q121Conf();
[two_sigmaMCMC, conf95MCMC] = q122Conf();

for p=1:NR_PARAMS
    h = figure
    plot(twoSigmaLap(p,:), [1,1],'m-x');
    hold on
    plot(two_sigmaParBoot(p,:), [2,2],'b-x');
    hold on
    plot(two_sigmaMCMC(p,:), [3,3],'g-x');
    hold on
    plot(conf95ParBoot(p,:), [4,4],'r-o');
    hold on
    plot(conf95MCMC(p,:), [5,5],'k-o');
    ylim ([0 6]);
    legend('2sigma Laplace','2 sigma Par Bootstrap', '2sigma MCMC', 'conf95ParBoot', 'conf95MCMC', 'location', 'northoutside');
    %legend('2 sigma Par Bootstrap', '2sigma MCMC', 'conf95ParBoot', 'conf95MCMC', 'location', 'northoutside');
    filename = sprintf('report/figures/q2/q123-p%d.eps', p);
    hgexport(h, filename);
end

%  diag cov - 2sigma range
%   -0.166775137824368
%   -0.000000000000000
%   -0.000000000016855
%   -0.000000000018668
%   -0.000000000025459

end

function twoSigma = LaplaceUncert(Avox, qhat, bvals, nr_iterations, startx, sigmaNoise)

[parameter_hat, ~, ~, Hessian] = fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx);

% correct by dividing with -2*sigma^2
cov = -inv(Hessian/(-2*sigmaNoise^2));

sigma = sqrt(diag(cov));

twoSigma = [(parameter_hat' -sigma), (parameter_hat' + sigma)];
twoSigma = twoSigma(1:3,:);

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
