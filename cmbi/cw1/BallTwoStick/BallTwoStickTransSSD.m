function [sumRes, S] = BallTwoStickTransSSD(x, Avox, bvals, qhat)

x = BallTwoStickTrans(x);
[sumRes, S] = BallTwoStickSSD(x, Avox, bvals, qhat);

end