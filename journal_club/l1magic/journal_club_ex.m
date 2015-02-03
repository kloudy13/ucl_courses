% l1qc_example.m
%
% Test out l1qc code (l1 minimization with quadratic constraint).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

% put optimization code in path if not already there
path(path, './Optimization');

% signal length
N = 1024;

% number of spikes to put down
T = 50;

% number of observations to make
K = 300;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));

% measurement matrix
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');
	
% noisy observations
sigma = 0.01;
e = sigma*randn(K,1);
y = A*x + e;

% initial guess = min energy
x0 = A'*y;

% take epsilon a little bigger than sigma*sqrt(K)
%epsilon =  sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K));
epsilon = 0.19 

tic
xp = l1qc_logbarrier(x0, A, [], y, epsilon, 0.1);
toc
norm1 = norm(xp - x)

[U, S, V] = svd(A);

perm = randperm(K);
Snew = zeros(size(S));
singVect = exp(-(1:K)/100);
singVect = singVect(perm);
Snew(1:K, 1:K) = diag(singVect);

Anew = U * Snew * V';

x0new = Anew'*y;

tic
xpnew = l1qc_logbarrier(x0new, Anew, [], y, epsilon, 0.1);
toc
norm2 = norm(xpnew - x)

yyy = 0;
% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1qc_logbarrier(x0, Afun, Atfun, y, epsilon, 1e-3, 50, 1e-8, 500);
% toc




