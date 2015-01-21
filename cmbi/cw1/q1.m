function q1()

format long 
NR_IMAGES = 33;
W = 112;
H = 112;

[dwis, qhat, bvals] = preprocessing(NR_IMAGES, W,H);

%q111(dwis, qhat, bvals);
%q112(dwis, qhat, bvals);
%q113(dwis, qhat, bvals);
q114(dwis, qhat, bvals);
%q116(dwis, qhat, bvals);
end

function [dwis, qhat, bvals] = preprocessing(NR_IMAGES, W,H)
fid = fopen('dwi.Bfloat','r','b');
dwis = fread(fid, 'float');
fclose(fid);

dwis = reshape(dwis, NR_IMAGES, W, H, 50);

qhat = load('grad_dirs.txt')';
bvals = 1000* sum(qhat .* qhat);
end

function [h] = eyeball(Avox, parameter_hat, bvals, qhat)

h = figure;
% Plot the actual data points
plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Predict measurements from model
predicted = BallStick(parameter_hat, bvals, qhat);
% Add the predictions to the plot
plot(predicted, ' rx', 'MarkerSize', 16, 'LineWidth', 4)
% Add labels and legend.
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');

end

function q111(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSD', startx, h, Avox, bvals, qhat)

eyeball(Avox, parameter_hat, bvals, qhat);

end

function q112(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% apply the inverse tranformations: sqrt and tangent
startx = [sqrt(7.5e+05) sqrt(3e-03) q1TransInv(2.5e-01) 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');


% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSDq112', startx, h, Avox, bvals, qhat);
RESNOM

% apply the transformations
[S0 d f theta phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
parameter_hat = [ S0^2 d^2 q1Trans(f) theta phi]

h = eyeball(Avox, parameter_hat, bvals, qhat);


end

function q113(dwis, qhat, bvals)

%Avox = dwis(:,52,62,25); %given voxel
Avox = dwis(:,23,40,18);
%Avox = dwis(:,70,64,14);
%Avox = dwis(:,111,111,14);


nr_iterations = 100;

[parameter_hat, minSSD] = fitVoxGlob1(Avox, qhat, bvals, nr_iterations)

h = eyeball(Avox, parameter_hat, bvals, qhat);

end

function [parameter_hat, minSSD] = fitVoxGlob1(Avox, qhat, bvals, nr_iterations)
% for use in q111-q114, uses a predefined starting point

% apply the inverse tranformations: sqrt and tangent
startx = [1.1e+05 2e-03 0.5 0 0];
startx = [sqrt(startx(1)) sqrt(startx(2)) q1TransInv(startx(3)) startx(4) startx(5)];

minSSD = inf;
minCounter = 0;
globTol = 0.1; % 0.1 recommended
bigSSDCount = 0;
minParHat = zeros(1,5);
sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
for i=1:nr_iterations
    
    sigma = eye(5);
    sigma(1,1) = sigmaScale*sqrt(7.5e+05);
    sigma(2,2) = sigmaScale*sqrt(3e-03);
    sigma(3,3) = sigmaScale*5;
    sigma(4,4) = sigAngleScale *pi;
    sigma(5,5) = sigAngleScale *pi;

    deltaX = mvnrnd(zeros(1,5),sigma);
    newStartX = startx + deltaX;
    
    % Define various options for the non-linear fitting algorithm
    h = optimset('MaxFunEvals', 20000, 'Algorithm', 'levenberg-marquardt',...
        'TolX', 1e-10, 'TolFun', 1e-10);


    % Now run the fitting
    % RESNOM is the value of the function at the solution found (parameter_hat)
    [parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fminunc('BallStickSSDq112', newStartX, h, Avox, bvals, qhat);
    %RESNOM

       
    parameter_hat(4:5) = mod(parameter_hat(4:5),2*pi);
    if (abs(minSSD - RESNOM) > globTol)
        abs(minSSD - RESNOM)
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
[S0 d f theta phi] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
parameter_hat = [ S0^2 d^2 q1Trans(f) theta phi];

end

function q114(dwis, qhat, bvals)

SLICE_NR=25;

[NR_IMAGES,W,H,~] = size(dwis);

%W = 5;
%H = 5;

mapS0 = zeros(W,H);
mapD = zeros(W,H);
mapF = zeros(W,H);
mapRESNORM = zeros(W,H);
mapTheta = zeros(W,H);
mapPhi = zeros(W,H);

nr_iterations = 2; % 2 recommended p(global_min) = 0.87

for w=1:W
    w
    for h=1:H   
        vox = dwis(:,w,h,SLICE_NR);
        [parameter_hat, mapRESNORM(w,h)] = fitVoxGlob1(vox, qhat, bvals, nr_iterations);
        [mapS0(w,h), mapD(w,h), mapF(w,h), mapTheta(w,h), mapPhi(w,h)] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
        
    end
end

save('q114', 'mapS0', 'mapD', 'mapF', 'mapRESNORM', 'mapTheta', 'mapPhi');

end

function q115(dwis, qhat, bvals)



end

function q116(dwis, qhat, bvals)

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit
startx = [7.5e+05 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm
h = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');
% S0 d f theta phi
lb = [0  , 0  , 0, -inf, -inf];
ub = [inf, inf, 1,  inf,  inf]; 
% Now run the fitting
% RESNOM is the value of the function at the solution found (parameter_hat)
[parameter_hat, RESNOM, EXITFLAG, OUTPUT] = fmincon('BallStickSSD', startx, [],[],[],[],lb, ub, [], h, Avox, bvals, qhat)


end

function y = q1Trans(x)
    y = 1/(1+exp(-x));
end

function x = q1TransInv(y)
    x = -log(1/y - 1)
end