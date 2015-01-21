function [h] = eyeball(Avox, parameter_hat, bvals, qhat)

h = figure;
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