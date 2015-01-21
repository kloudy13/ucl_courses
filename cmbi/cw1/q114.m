function q114(dwis, qhat, bvals)

[dwis, qhat, bvals] = q1Preprocessing();

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
startx = [1.1e+05 2e-03 0.5 0 0];

for w=1:W % errors on w=99 h=63
    w
    for h=1:H   
        vox = dwis(:,w,h,SLICE_NR);
        succeeded = false;
        while ~succeeded
          try

            [parameter_hat, mapRESNORM(w,h)] = fitVoxGlob1(vox, qhat, bvals, nr_iterations, startx);
            succeeded = true;
          catch
          end
        end
        [mapS0(w,h), mapD(w,h), mapF(w,h), mapTheta(w,h), mapPhi(w,h)] = deal(parameter_hat(1),parameter_hat(2),parameter_hat(3),parameter_hat(4),parameter_hat(5));
        
    end
end

save('q114', 'mapS0', 'mapD', 'mapF', 'mapRESNORM', 'mapTheta', 'mapPhi');

end