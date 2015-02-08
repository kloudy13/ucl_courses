function y = NoddiSimpleTransInv(x)
assert(x(1) >= 0 && x(2) >= 0);
y(1) = sqrt(x(1)); %S0
y(2) = sqrt(x(2)); % d
y(3) = q1TransInv(x(3)); % f
y(4) = x(4); % theta
y(5) = x(5); % phi
y(6) = q1TransInv(x(6)); % k

end