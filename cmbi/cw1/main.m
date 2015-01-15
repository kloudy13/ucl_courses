function q1()
end

function q111()

fid = fopen('dwi.Bfloat','r','b');
dwis = fread(fid, 'float');
fclose(fid);

NR_IMAGES = 33;
W = 112;
H = 112;

dwis = reshape(dwis, NR_IMAGES, W, H, 50);

%imshow(squeeze(dwis(1,:,:,25)), []);
%imshow(squeeze(dwis(2,:,:,25)), []);

qhat = load('grad_dirs.txt')';
bvals = 1000* sum(qhat .* qhat);

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
%startx = [7.5e+05 3e-03 2.5e-01 0 0]; % correct one
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSD', startx, h, Avox, bvals, qhat)


%q1.1.1

figure;
% Plot the actual data points
plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Predict measurements from model
predicted = BallStick(parameter_hat, bvals, qhat);
% Add the predictions to the plot
plot(predicted, ' rx', 'MarkerSize', 16, 'LineWidth', 4)
% Add labels and legend.
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');

end

