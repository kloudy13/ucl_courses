function q113(dwis, qhat, bvals)
[dwis, qhat, bvals] = q1Preprocessing();

%Avox = dwis(:,52,62,25); %given voxel
Avox = dwis(:,63,40,18);
%Avox = dwis(:,70,64,14);
%Avox = dwis(:,111,111,14);
%Avox = dwis(:,99,63,25);
%Avox = dwis(:,32,15,25);

nr_iterations = 100;

startx = [1.1e+05 2e-03 0.5 0 0];
[parameter_hat, minSSD] = fitVoxGlob1(Avox, qhat, bvals, nr_iterations, startx)

h = eyeball(Avox, parameter_hat, bvals, qhat);

end