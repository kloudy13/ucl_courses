function [sumRes, S] = NoddiSimpleSSD(x, Avox, bvals, qhat)

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
k = x(6);
%lam1 = x(7);
%lam2 = x(8);

% Synthesize the signals
mu = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

NR_MEASUREMENTS = length(qhat);
NR_SAMPLES = 2000;
SIsamples = zeros(NR_SAMPLES, NR_MEASUREMENTS);
n_samples = RandSampleSphere(NR_SAMPLES,'uniform');
for sample_nr=1:NR_SAMPLES
  n = n_samples(sample_nr,:);
  qDotNSquared = (sum(qhat .* repmat(n, [NR_MEASUREMENTS 1])')).^2;
  SIsamples(sample_nr,:) = watson(k,n,mu)*exp(-bvals * diff .*  qDotNSquared); % intra-cellular diffusion
end

n_samplesM = repmat(n_samples, [1 1 NR_MEASUREMENTS]);
qhatS = repmat(qhat, [1 1 NR_SAMPLES]);
n_samplesM = permute(n_samplesM,[2 3 1]);
qDotNSquared2 = squeeze(sum(qhatS .* n_samplesM, 1).^2);
bvalsS = repmat(bvals, [NR_SAMPLES 1])';
watsonM = repmat(watson(k,n_samples,mu), [NR_MEASUREMENTS 1])';
SIsamples2 = *exp(-bvals * diff .*  qDotNSquared);

Si = mean(SIsamples); % multiply by the surface of the sphere

%Se = exp(-bvals.*(lam2 + (lam1 - lam2)*muDotNSquared)); % extra-cellular diffusion
Se = exp(-bvals*diff); % extra-cellular diffusion

S = S0*(f*Si +(1-f)*Se);

% Compute the sum of squared differences
sumRes = sum((Avox - S').^2);

end



