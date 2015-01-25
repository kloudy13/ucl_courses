function q122()

[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25);

NR_PARAMS = 3;
NR_SAMPLES = 200;
sigmaQ = [10^3, 0.0001 0.1];

params = [1.1e+05 2e-03 0.5 0 0];

samples = zeros(NR_PARAMS,NR_SAMPLES);
for p=1:NR_PARAMS
  [samples(p,:), acc_rate] = MCMC(params, p, Avox, qhat, bvals, NR_SAMPLES, sigmaQ);
  acc_rate
end

for p=1:NR_PARAMS
    h = q12calcUncertainty(samples(p,:));  
    filename = sprintf('report/figures/q2/q122-p%d.eps', p);
    hgexport(h, filename);
end

end

function [samples, acc_rate] = MCMC(params, parIndex, Avox, qhat, bvals, NR_SAMPLES, sigmaQ)

samples = zeros(NR_SAMPLES, 1);

acceptanceCount = 0;

xOld = params(parIndex);
lik_old = BallStickSSD(params, Avox, bvals, qhat);
for t=1:NR_SAMPLES
  xNew = normrnd(xOld, sigmaQ(parIndex));
  paramsNew = params;
  paramsNew(parIndex) = xNew;
  a = BallStickSSD(paramsNew, Avox, bvals, qhat)/ lik_old;
  if (a > rand)
    samples(t) = xNew;
    lik_old = BallStickSSD(paramsNew, Avox, bvals, qhat);
    acceptanceCount = acceptanceCount + 1;
  else
    samples(t) = xOld;
  end

  
end

acc_rate = acceptanceCount / NR_SAMPLES;

end