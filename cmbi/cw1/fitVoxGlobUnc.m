function [parameter_hat, minSSD, minCounter] = fitVoxGlobUnc(Avox, qhat, bvals, nr_iterations, startx, sigma, fminuncOptions, globTol, model)
% for use in q111-q114, uses a predefined starting point

% apply the inverse tranformations: sqrt and sigmoid
startxTrans = startx;
startxTrans(1:3) = [sqrt(startx(1)) sqrt(startx(2)) q1TransInv(startx(3))];


minCounter = 0;
bigSSDCount = 0;

[minParHat, minSSD] = fminunc(model, startxTrans, fminuncOptions, Avox, bvals, qhat);

NR_PARAMS = length(startx);

%newXTrans = newStartX;
%newXTrans(1:3) = [sqrt(newStartX(1)) sqrt(newStartX(2)) q1TransInv(newStartX(3))];

for i=1:nr_iterations
    i
    deltaX = mvnrnd(zeros(1,NR_PARAMS),sigma);
    newXTrans = startxTrans + deltaX;
    

    succeeded = false;
    while ~succeeded
      try
        % Now run the fitting ... if the gradient method fails because of
        % approximation errors on the Hessian, try again. Happends very rarely
        [parameter_hat, RESNOM] = fminunc(model, newXTrans, fminuncOptions, Avox, bvals, qhat);
        succeeded = true;
      catch
        newXTrans
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
    end
    

    
end
parameter_hat = minParHat;
% apply the transformations
[S0, d, f, theta, phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
parameter_hat(1:3) = [ S0^2 d^2 q1Trans(f)];
%minSSD
%minCounter
end