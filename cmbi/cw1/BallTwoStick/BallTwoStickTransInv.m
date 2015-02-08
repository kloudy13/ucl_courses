function y = BallTwoStickTransInv(x)
assert(x(1) >= 0 && x(2) >= 0);
y(1) = sqrt(x(1)); %S0
y(2) = sqrt(x(2)); %d
y(3) = sqrt(x(3)); %f1
y(4) = sqrt(x(4)); %f2
y(5) = x(5); %theta1
y(6) = x(6); %phi1
y(7) = x(7); %theta2
y(8) = x(8); %phi2


end