function q114(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

SLICE_NR=25;

[NR_IMAGES,W,H,~] = size(dwis);

mapS0 = zeros(W,H);
mapD = zeros(W,H);
mapF = zeros(W,H);
mapRESNORM = zeros(W,H);
mapTheta = zeros(W,H);
mapPhi = zeros(W,H);

nr_iterations = 2; % 2 recommended p(global_min) = 0.87
startx = [1.1e+05 2e-03 0.5 0 0];
fminuncOptions = optimoptions(@fminunc,'Algorithm','quasi-newton', 'MaxFunEvals', 20000,'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off'); 

sigAngleScale = 2; % 2 recommended
sigmaScale = 10; % 10 recommended
globTol = 0.1; % 0.1 recommended

sigma = eye(5);
sigma(1,1) = sigmaScale*sqrt(7.5e+05);
sigma(2,2) = sigmaScale*sqrt(3e-03);
sigma(3,3) = sigmaScale*5;
sigma(4,4) = sigAngleScale *pi;
sigma(5,5) = sigAngleScale *pi;
model = 'BallStickSSDq112';

for w=1:W % errors on w=99 h=63
    w
    for h=1:H   
        vox = dwis(:,w,h,SLICE_NR);
        [parameter_hat, mapRESNORM(w,h)] = fitVoxGlobUnc(vox, qhat, bvals, nr_iterations, startx, sigma, fminuncOptions, globTol, model);
        [mapS0(w,h), mapD(w,h), mapF(w,h), mapTheta(w,h), mapPhi(w,h)] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
        
    end
end

save('q114', 'mapS0', 'mapD', 'mapF', 'mapRESNORM', 'mapTheta', 'mapPhi');

end