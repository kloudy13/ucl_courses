function [sumRes, S] = ZeppTwoStickTransSSD(x, Avox, bvals, qhat)

x = ZeppTwoStickTrans(x);
[sumRes, S] = ZeppTwoStickSSD(x, Avox, bvals, qhat);

end