function [sumRes, S] = q13ZeppelinStickTransSSD(x, Avox, bvals, qhat)

x = q3ZeppTrans(x);
[sumRes, S] = q13ZeppelinStickSSD(x, Avox, bvals, qhat);

end