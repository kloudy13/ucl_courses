function q142()

[signals, bvals, qhat] = q13preprocessing();

load('q131BallStick.mat');
sigma = 0.04;

NR_SHELLS = 24;
NR_MEASUREMENTS=48;

Fshells = zeros(NR_SHELLS,3,3); 

for shell_nr=0:NR_SHELLS-1
  indices = shell_nr*NR_MEASUREMENTS+1:(shell_nr+1)*NR_MEASUREMENTS;
  %sigmalsShell = sigmals(indices);
  F = Fisher(paramsBallStick, signals(indices), bvals(:,indices), qhat(:,indices), sigma);
  Fshells(shell_nr+1,:,:) = F(1:3, 1:3);
end

[bestScoreA, bestShellA, scoresA] = Aoptimality(Fshells);
[bestScoreT, bestShellT, scoresT] = Toptimality(Fshells);

end

function [bestScore, bestShell, scores] = Aoptimality(Fshells)

NR_SHELLS = size(Fshells, 1);
scores = zeros(NR_SHELLS,1);
for shell_nr = 1:NR_SHELLS
  scores(shell_nr) = trace(inv(squeeze(Fshells(shell_nr,:,:))));
end

[bestScore, bestShell] = min(scores);

end

function [bestScore, bestShell, scores] = Toptimality(Fshells)

NR_SHELLS = size(Fshells, 1);
scores = zeros(NR_SHELLS,1);
for shell_nr = 1:NR_SHELLS
  scores(shell_nr) = trace(squeeze(Fshells(shell_nr,:,:)));
end

[bestScore, bestShell] = max(scores);

end

function F = Fisher(paramsBallStick, signals, bvals, qhat, sigma)

dSdS0 = BallStickDerivSigS0(paramsBallStick, bvals, qhat);
dSdD = BallStickDerivSigD(paramsBallStick, bvals, qhat);
dSdF = BallStickDerivSigF(paramsBallStick, bvals, qhat);
dSdTheta = BallStickDerivSigTheta(paramsBallStick, bvals, qhat);
dSdPhi = BallStickDerivSigPhi(paramsBallStick, bvals, qhat);

dS = [dSdS0; dSdD; dSdF; dSdTheta; dSdPhi]';

NR_IMAGES = length(signals);
NR_PARAMS = length(paramsBallStick);
F = zeros(NR_PARAMS, NR_PARAMS);
for k=1:NR_IMAGES
  F = F + dS(k,:)' * dS(k,:);
end

F = F / sigma^2;
F = F ./ (paramsBallStick' * paramsBallStick);

end