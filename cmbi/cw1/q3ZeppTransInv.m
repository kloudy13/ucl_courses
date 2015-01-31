function [S0, diff, f, theta, phi, lam1, lam2] = q3ZeppTransInv(x)
S0 = sqrt(x(1));
diff = sqrt(x(1));
f = q1TransInv(x(3));
theta = x(4);
phi = x(5);
lam1 = x(6);
lam2 = x(7);

end