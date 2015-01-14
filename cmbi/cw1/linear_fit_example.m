% Gradient directions in qhat (0 for b=0 measurements).
% Load data into A (see coursework)
A = dwis(:,66,58,25);
% Construct the design matrix
b = 1000;
G = [ones(1, length(qhat)); -b*sum(qhat.*qhat)]';
% Estimate parameters
x = pinv(G)*log(A); % Or, maybe better: x = G\log(A);
logS0 = x(1);
adc = x(2);


figure;
% Plot the actual data points
plot(A, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Predict measurements from model
model = exp(G*x);
% Add the predictions to the plot
plot(model, ' rx', 'MarkerSize', 16, 'LineWidth', 4)
% Add labels and legend.
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');