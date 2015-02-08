function y = ZeppTwoStickTrans(x)
y(1) = x(1)^2; %S0
y(2) = x(2)^2; % d
y(3) = x(3)^2/(x(3)^2+x(4)^2+0.33); % f1
y(4) = x(4)^2/(x(3)^2+x(4)^2+0.33); % f2
y(5) = x(5); % theta1
y(6) = x(6); % phi1
y(7) = x(7); % theta2
y(8) = x(8); % phi2
y(9) = x(9)^2; %lam1
y(10) = x(10)^2; %lam2


end