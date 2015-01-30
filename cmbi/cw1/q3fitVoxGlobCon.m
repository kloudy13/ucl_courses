function [parameter_hat, minSSD, minCounter, globMinHessian] = q3fitVoxGlobCon(Avox, qhat, bvals, nr_iterations, startx, lb, ub, sigma, fminconOptions, model)
% for use in q111-q114, uses a predefined starting point

% apply the inverse tranformations: sqrt and tangent
    
% Define various options for the non-linear fitting algorithm
%options = optimset('MaxFunEvals', 20000, 'Algorithm', 'quasi-newton',...
%    'TolX', 1e-10, 'TolFun', 1e-10);

%options = optimset('MaxFunEvals', 20000, 'Algorithm', 'active-set');


[~, minSSD] = fmincon(model, startx, [],[],[],[],lb, ub, [], fminconOptions, Avox, bvals, qhat);

minCounter = 0;
globTol = 0.1; % 0.1 recommended
bigSSDCount = 0;
minParHat = startx;
NR_PARAMS = length(startx);

for i=1:nr_iterations

    deltaX = mvnrnd(zeros(1,NR_PARAMS),sigma.^2);
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
        [parameter_hat, RESNOM, exitflag, output, ~, ~, Hessian] = fmincon(model, newStartX, [],[],[],[],lb, ub, [], fminconOptions, Avox, bvals, qhat);
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