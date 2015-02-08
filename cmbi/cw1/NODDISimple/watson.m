function x = watson(k, n_samples, mu)
% n is a nx3 vector of possible fibre directions over which we integrate by
% sampling
% mu is a 3x1 vector of the mean direction of the Watson distribution
% k is the concentration parameter

NR_SAMPLES = size(n_samples,1);
muS = repmat(mu, [NR_SAMPLES 1]);
x = KummerComplex(1/2, 3/2, k)^(-1) * exp( k * sum(n_samples .* muS, 2)^2);

end