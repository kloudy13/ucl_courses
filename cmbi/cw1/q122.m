function [two_sigma, conf95] = q122()

[dwis, qhat, bvals] = q1Preprocessing();

%Avox = dwis(:,52,62,25);
%Avox = dwis(:,63,40,18);
Avox = dwis(:,70,64,14);

NR_PARAMS = 3;
NR_SAMPLES = 100000;
NR_BURN_IN = NR_SAMPLES / 10;
sigmaQ = [10^4, 0.0001 0.05];
sigmaNoise = 5000;

params = [1.132e+05 1.534e-03 0.575 2.106 6.170];
%params = 10^5 * [1.132129187022450   0.000000015341202   0.000005754246427   0.000021065381746   0.000061701371023];
%params = [1.1e+05 7.534e-03 0.175 0 0];

% samples = zeros(NR_PARAMS,NR_SAMPLES+NR_BURN_IN);
% for p=1:NR_PARAMS
%   [samples(p,:), acc_rate] = MCMC(params, p, Avox, qhat, bvals, NR_SAMPLES+NR_BURN_IN, sigmaQ, sigmaNoise);
%   mean(samples(p,:))
%   plot(samples(p,:))
%   acc_rate
% end

covQ = eye(5);
covQ(1,1) = 10^4;
covQ(2,2) = 0.00005;
covQ(3,3) = 0.01;
covQ(4,4) = 0.1;
covQ(5,5) = 0.1;
% raise to power 2 since mvnrnd requires covariance and not std deviations
covQ = covQ .^2;
sigmaNoise = 10000;


[samples, acc_rate] = MCMCall3(params, Avox, qhat, bvals, NR_SAMPLES+NR_BURN_IN, covQ, sigmaNoise);

acc_rate

for p=1:NR_PARAMS
    plot(samples(p, NR_BURN_IN:end));
    mean(samples(p, NR_BURN_IN:end))
end

samples = samples(:, NR_BURN_IN:end);
save('q122-vox3.mat', 'samples');

two_sigma = zeros(NR_PARAMS, 2);
conf95 = zeros(NR_PARAMS, 2);
for p=1:NR_PARAMS
    [h, two_sigma(p,:), conf95(p, :)] = q12calcUncertainty(samples(p,:), 1);  
    filename = sprintf('report/figures/q2/q122-p%d.eps', p);
    %hgexport(h, filename);
end

end


function [samples, acc_rate] = MCMCall3(params, Avox, qhat, bvals, NR_SAMPLES, covQ, sigmaNoise)

samples = zeros(3, NR_SAMPLES);

acceptanceCount = 0;

paramsOld = params;
ssd_old = BallStickSSD(params, Avox, bvals, qhat);
paramsNew = params;
for t=1:NR_SAMPLES
  paramsNew = mvnrnd(paramsOld, covQ); % sample new data point
  ssd_new = BallStickSSD(paramsNew, Avox, bvals, qhat); % compute new ssd
  a = (ssd_old - ssd_new)/(2*sigmaNoise^2); % calc acceptance ratio
  if (a > log(rand))
    % accept the sample
    samples(:,t) = paramsNew(1:3);
    ssd_old = ssd_new;
    acceptanceCount = acceptanceCount + 1;
    paramsOld = paramsNew;
  else
    % reject the sample
    samples(:,t) = paramsOld(1:3);
  end

  
end

acc_rate = acceptanceCount / NR_SAMPLES;

end


function [samples, acc_rate] = MCMC(params, parIndex, Avox, qhat, bvals, NR_SAMPLES, sigmaQ, sigmaNoise)

samples = zeros(NR_SAMPLES, 1);

acceptanceCount = 0;

xOld = params(parIndex);
ssd_old = BallStickSSD(params, Avox, bvals, qhat);
paramsNew = params;
for t=1:NR_SAMPLES
  xNew = normrnd(xOld, sigmaQ(parIndex)); % sample new data point
  paramsNew(parIndex) = xNew; % update paramsNew vector with new point
  ssd_new = BallStickSSD(paramsNew, Avox, bvals, qhat); % compute new ssd
  a = (ssd_old - ssd_new)/(2*sigmaNoise^2); % calc acceptance ratio
  if (a > log(rand))
    % accept the sample
    samples(t) = xNew;
    ssd_old = ssd_new;
    acceptanceCount = acceptanceCount + 1;
    xOld = xNew;
  else
    % reject the sample
    samples(t) = xOld;
  end

  
end

acc_rate = acceptanceCount / NR_SAMPLES;

end
