function [sumRes, S] = NoddiSimpleTransSSD(x, Avox, bvals, qhat)

x = NoddiSimpleTrans(x);
[sumRes, S] = NoddiSimpleSSD(x, Avox, bvals, qhat);

end