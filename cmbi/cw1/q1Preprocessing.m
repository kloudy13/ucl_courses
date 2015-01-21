function [dwis, qhat, bvals] = q1Preprocessing()

format long 
NR_IMAGES = 33;
W = 112;
H = 112;

fid = fopen('dwi.Bfloat','r','b');
dwis = fread(fid, 'float');
fclose(fid);

dwis = reshape(dwis, NR_IMAGES, W, H, 50);

qhat = load('grad_dirs.txt')';
bvals = 1000* sum(qhat .* qhat);
end