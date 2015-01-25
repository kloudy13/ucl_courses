function h = q12calcUncertainty(parameter_sample)
% input is a 1D parameter vector, respresenting samples for one of the 5
% parameters

mx = mean(parameter_sample);
sx = std(parameter_sample);
N = length(parameter_sample);

% Approximative
upper_limit=mx+1.96*sx/sqrt(N);
lower_limit=mx-1.96*sx/sqrt(N);
%plot(1:N,parameter_sample,'o',1:N,upper_limit,1:N,lower_limit)

% 0.025-0.975
sorted_samples = sort(parameter_sample);
[lower_limit, upper_limit] = deal(sorted_samples(floor(0.025*N)),sorted_samples(ceil(0.975*N)));  

[lower_limit, upper_limit]
h = figure
hist(parameter_sample, 20)
maxY = ylim;
hold on
plot([lower_limit, upper_limit], 1.05*[maxY(2), maxY(2)],'k--o')
%plot([lower_limit, upper_limit], [70,70],'k--o')
hold on
plot([mx - sx, mx + sx], 1.1*[maxY(2), maxY(2)], 'r--o')
legend('p(x|A)','95% range','2 sigma range');
end