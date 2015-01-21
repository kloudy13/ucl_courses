function q115(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

%voxSig = dwis(:,52,62,25);
voxSig = dwis(:,63,40,18);
%voxSig = dwis(:,70,64,14);
%voxSig = dwis(:,111,111,14);
%voxSig = dwis(:,99,63,25);
%voxSig = dwis(:,32,15,25);

[parameter_hat, minSSD] = q1FitVoxLin(voxSig, qhat, bvals)

end