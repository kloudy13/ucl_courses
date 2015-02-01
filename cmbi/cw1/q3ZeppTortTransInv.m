function y = q3ZeppTortTransInv(x)
assert(x(1) >= 0 && x(2) >= 0 && x(6) >=0);
y(1) = sqrt(x(1));
y(2) = sqrt(x(2));
y(3) = q1TransInv(x(3));
y(4) = x(4);
y(5) = x(5);
y(6) = sqrt(x(6));

end