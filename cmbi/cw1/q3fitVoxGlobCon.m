function [parameter_hat, minSSD, minCounter, globMinHessian] = q3fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx)
% for use in q111-q114, uses a predefined starting point

% apply the inverse tranformations: sqrt and tangent
    
% Define various options for the non-linear fitting algorithm
%options = optimset('MaxFunEvals', 20000, 'Algorithm', 'quasi-newton',...
%    'TolX', 1e-10, 'TolFun', 1e-10);

%options = optimset('MaxFunEvals', 20000, 'Algorithm', 'active-set');


lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 
options = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');

[~, minSSD] = fmincon('BallStickSSD', startx, [],[],[],[],lb, ub, [], options, Avox, bvals, qhat);

minCounter = 0;
globTol = 0.1; % 0.1 recommended
bigSSDCount = 0;
minParHat = startx;
sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended

for i=1:nr_iterations
    
    sigma = eye(5);
    sigma(1,1) = 0.15;
    sigma(2,2) = 1e-09;
    sigma(3,3) = 0.15;
    sigma(4,4) = pi;
    sigma(5,5) = pi;

    deltaX = mvnrnd(zeros(1,5),sigma.^2);
    newStartX = startx + deltaX;
    startingSSD = BallStickSSD(newStartX, Avox, bvals, qhat)

    if(newStartX(1) < 0 || newStartX(2) < 0 || abs(newStartX(3) - 0.5) > 0.5)
       display('newStartX out of bounds'); 
    end
    
    succeeded = false;
    while ~succeeded
      try
        % Now run the fitting ... if the gradient method fails because of
        % approximation errors on the Hessian, try again. Happends very rarely
        [parameter_hat, RESNOM, exitflag, output, ~, ~, Hessian] = fmincon('BallStickSSD', newStartX, [],[],[],[],lb, ub, [], options, Avox, bvals, qhat);
        succeeded = true;
      catch
      end
    end
    
    
    parameter_hat(4:5) = mod(parameter_hat(4:5),2*pi);
    if (abs(minSSD - RESNOM) > globTol)
        %abs(minSSD - RESNOM)
        bigSSDCount = bigSSDCount+1;
        sum(abs(parameter_hat - minParHat))
    end
    
    if(abs(minSSD - RESNOM) <= globTol)
       % found a previous minimum
       %minSSD
       minCounter = minCounter + 1; 
    end
    
    if (abs(minSSD - RESNOM) > globTol && minSSD > RESNOM)
        % found a new global minimum
       minSSD = RESNOM;
       %minSSD
       minCounter = 0;
       minParHat = parameter_hat;
       globMinHessian = Hessian;
    end
    

    
end
parameter_hat = minParHat;
[S0, d, f, theta, phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
%minSSD
%minCounter
end