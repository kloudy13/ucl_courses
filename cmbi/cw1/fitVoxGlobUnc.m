function [parameter_hat, minSSD, minCounter] = fitVoxGlobUnc(Avox, qhat, bvals, nr_iterations, startx, sigma, fminuncOptions, globTol, model)
% for use in q111-q114, uses a predefined starting point

% apply the inverse tranformations: sqrt and sigmoid
startxTrans = startx;
startxTrans(1:3) = [sqrt(startx(1)) sqrt(startx(2)) q1TransInv(startx(3))];


minCounter = 0;
bigSSDCount = 0;

[minParHatTrans, minSSD] = fminunc(model, startxTrans, fminuncOptions, Avox, bvals, qhat);

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
        [parameter_hatTrans, RESNOM] = fminunc(model, newXTrans, fminuncOptions, Avox, bvals, qhat);
        succeeded = true;
      catch
        [newXTrans(1)^2,newXTrans(2)^2, q1Trans(newXTrans(3))] 
      end
    end
    
    
    parameter_hatTrans(4:5) = mod(parameter_hatTrans(4:5),2*pi);
    if (abs(minSSD - RESNOM) > globTol)
        %abs(minSSD - RESNOM)
        bigSSDCount = bigSSDCount+1;
        sum(abs(parameter_hatTrans - minParHatTrans))
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
       minParHatTrans = parameter_hatTrans;
    end
    

    
end
parameter_hat = minParHatTrans;
% apply the transformations
[S0, d, f, theta, phi] = deal(parameter_hatTrans(1),parameter_hatTrans(2),parameter_hatTrans(3),parameter_hatTrans(4),parameter_hatTrans(5));
parameter_hat(1:3) = [ S0^2 d^2 q1Trans(f)];
%minSSD
%minCounter
end