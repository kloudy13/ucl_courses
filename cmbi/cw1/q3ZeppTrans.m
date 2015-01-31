function [S0, diff, f, theta, phi, lam1, lam2] = q3ZeppTrans(x)
S0 = x(1)^2;
diff = x(2)^2;
f = q1Trans(x(3));
theta = x(4);
phi = x(5);
lam1 = x(6)^2;
lam2 = x(7)^2;

end