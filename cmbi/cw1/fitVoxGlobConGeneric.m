function [parameter_hat, minSSD, minCounter, globMinHessian] = fitVoxGlobConGeneric(Avox, qhat, bvals, nr_iterations, startx, lb, ub,sigma,globTol, fmincon_options,model)


minCounter = 0;

bigSSDCount = 0;
% sigAngleScale = 2; % 2 recommended
% sigmaScale = 10; % 10 recommended


[minParHat, minSSD, ~, ~, ~, ~, globMinHessian] = fmincon(model, startx, [],[],[],[],lb, ub, [], fmincon_options, Avox, bvals, qhat);


for i=1:nr_iterations


    deltaX = mvnrnd(zeros(1,5),sigma);
    newStartX = startx + deltaX;


    succeeded = false;
    while ~succeeded
      try
        % Now run the fitting ... if the gradient method fails because of
        % approximation errors on the Hessian, try again. Happends very rarely
        [parameter_hat, RESNOM, ~, ~, ~, ~, Hessian] = fmincon(model, newStartX, [],[],[],[],lb, ub, [], fmincon_options, Avox, bvals, qhat);
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
% [S0, d, f, theta, phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
%minSSD
%minCounter
end