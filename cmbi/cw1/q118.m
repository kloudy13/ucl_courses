function q118()
[dwis, qhat, bvals] = q1Preprocessing();

Avox = dwis(:,52,62,25); %given voxel
%Avox = dwis(:,63,40,18);
%Avox = dwis(:,70,64,14);
%Avox = dwis(:,111,111,14);
%Avox = dwis(:,99,63,25);
%Avox = dwis(:,32,15,25);

nr_iterations = 30;

startx = [1.1e+05 1.5e-03 0.5 0 0];
[parameter_hat, minSSD] = fitVoxGlobRician(Avox, qhat, bvals, nr_iterations, startx)

predicted = BallStick(parameter_hat, bvals, qhat);
h = eyeball(Avox, predicted, bvals, qhat);

end