function [two_sigma, conf95] = q122()

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25)';

NR_PARAMS = 3;
NR_SAMPLES = 10000;
sigmaQ = [10^3, 0.001 0.1];
sigmaNoise = 10000;

params = [1.1e+05 1.5e-03 0.5 0 0];

samples = zeros(NR_PARAMS,NR_SAMPLES);
for p=1:NR_PARAMS
  [samples(p,:), acc_rate] = MCMC(params, p, Avox, qhat, bvals, NR_SAMPLES, sigmaQ, sigmaNoise);
  mean(samples(p,:))
  plot(samples(p,:))
  acc_rate
end


save('q122.mat', 'samples');

two_sigma = zeros(NR_PARAMS, 2);
conf95 = zeros(NR_PARAMS, 2);
for p=1:NR_PARAMS
    [h, two_sigma(p,:), conf95(p, :)] = q12calcUncertainty(samples(p,:), 1);  
    filename = sprintf('report/figures/q2/q122-p%d.eps', p);
    hgexport(h, filename);
end

end

function [samples, acc_rate] = MCMC(params, parIndex, Avox, qhat, bvals, NR_SAMPLES, sigmaQ, sigmaNoise)

samples = zeros(NR_SAMPLES, 1);

acceptanceCount = 0;

xOld = params(parIndex);
lik_old = BallStickProb(params, Avox, bvals, qhat, sigmaNoise);
for t=1:NR_SAMPLES
  xNew = normrnd(xOld, sigmaQ(parIndex));
  paramsNew = params;
  paramsNew(parIndex) = xNew;
  lik_new = BallStickProb(paramsNew, Avox, bvals, qhat, sigmaNoise);
  a = lik_new/ lik_old;
  if (a > rand)
    % accept the sample
    samples(t) = xNew;
    lik_old = lik_new;
    acceptanceCount = acceptanceCount + 1;
    xOld = xNew;
  else
    % reject the sample
    samples(t) = xOld;
  end

  
end

acc_rate = acceptanceCount / NR_SAMPLES;

end