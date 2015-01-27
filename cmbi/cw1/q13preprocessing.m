function [signals, bvals, qhat] = q13preprocessing()

datafile = 'challengeOPEN.txt';
fid = fopen(datafile, 'r', 'b');
% Read in the header
A = fgetl(fid);
	
% Read in the data
A = fscanf(fid, '%f', [8, inf]);
fclose(fid);
% Create the protocol
signals = A(1,:)';
grad_dirs = A(2:4,:);
G = A(5,:);
bigD = A(6,:);
smallD = A(7,:);
TE = A(8,:);
GAMMA = 2.675987E8;

bvals = ((GAMMA*smallD.*G).^2).*(bigD-smallD/3);

% set qhat as the gradient directions
qhat = grad_dirs;
end