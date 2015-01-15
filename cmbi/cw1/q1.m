function q1()

NR_IMAGES = 33;
W = 112;
H = 112;

[dwis, qhat, bvals] = preprocessing(NR_IMAGES, W,H);

%q111(dwis, qhat, bvals);
q112(dwis, qhat, bvals);
%q116(dwis, qhat, bvals);
end

function [dwis, qhat, bvals] = preprocessing(NR_IMAGES, W,H)
fid = fopen('dwi.Bfloat','r','b');
dwis = fread(fid, 'float');
fclose(fid);

dwis = reshape(dwis, NR_IMAGES, W, H, 50);

qhat = load('grad_dirs.txt')';
bvals = 1000* sum(qhat .* qhat);
end

function q111(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSD', startx, h, Avox, bvals, qhat)

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

function q112(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% apply the inverse tranformations: sqrt and tangent
startx = [sqrt(7.5e+05) sqrt(3e-03) tan(2.5e-01) 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSDq112', startx, h, Avox, bvals, qhat);
RESNOM

% apply the transformations
[S0 d f theta phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
parameter_hat = [ S0^2 d^2 atan(f) theta phi]

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

function q116(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi
lb = [0, 0, 0, -inf, -inf];
ub = [inf, inf, 1, inf, inf]; 
% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fmincon('BallStickSSD', startx, [],[],[],[],lb, ub, [], h, Avox, bvals, qhat)

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
