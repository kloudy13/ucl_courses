function [h] = q3eyeball(Avox, predicted, bvals, qhat)

h = figure('Position', [300, 0, 700, 500]);;
% Plot the actual data points
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
% Add the predictions to the plot
plot(predicted, ' rx', 'MarkerSize', 6, 'LineWidth', 2)
% Add labels and legend.
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 15);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model', 'Location','northoutside');

end