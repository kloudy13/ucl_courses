function [sumRes, S] = q13ZeppelinStickTortTransSSD(x, Avox, bvals, qhat)

x = q3ZeppTortTrans(x);
[sumRes, S] = q13ZeppelinStickTortSSD(x, Avox, bvals, qhat);

end