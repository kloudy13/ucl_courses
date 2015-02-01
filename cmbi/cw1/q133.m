function q133()

[signals, bvals, qhat] = q13preprocessing();

K = length(signals);
load('q131BallStick.mat');
aicBallStick =calcAIC(5, K, SSDDBallStick)
bicBallStick =calcBIC(5, K, SSDDBallStick)

% Zeppelin-Stick
load('q132ZeppStickUnc.mat');
aicZeppStick =calcAIC(7, K, SSDZeppStick)
bicZeppStick =calcBIC(7, K, SSDZeppStick)

% Zeppelin-Stick with tortuosity
load('q132ZeppStickTortUnc.mat');
aicZeppStickTort =calcAIC(6, K, SSDZeppStickTort)
bicZeppStickTort =calcBIC(6, K, SSDZeppStickTort)


end

function aic = calcAIC(N, K, SSD)
% N is the nr of parameters 
% K is the number of data points
aic = 2*N + K*log(SSD/K);
end

function bic = calcBIC(N, K, SSD)
% N is the nr of parameters 
% K is the number of data points
bic = N*log(K) + K*log(SSD/K);
end