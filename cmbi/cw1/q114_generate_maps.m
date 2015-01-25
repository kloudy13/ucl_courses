function q114_generate_maps()

load('q114.mat')

maskThreshMin = 5*10^8;
maskThreshMax = 4*10^9;
mask = (mapRESNORM > maskThreshMin) .* (mapRESNORM < maskThreshMax);
mapRESNORM = mapRESNORM .* mask;

figure('name','Fbef');
imshow(mapF, [min(mapF(:)), max(mapF(:))])

mapD = mapD .* mask;
mapF = mapF .* mask;
mapS0 = mapS0 .* mask;

figure('name','Fafter');
imshow(mapF, [min(mapF(:)), max(mapF(:))])

fibDirXMap = mapF .* cos(mapPhi) .* sin(mapTheta);
fibDirYMap = mapF .* sin(mapPhi) .* sin(mapTheta);
fibDirXMap = fibDirXMap .* mask;
fibDirYMap = fibDirYMap .* mask;

h = figure('name','S0');
imshow(mapS0, [min(mapS0(:)), max(mapS0(:))])
hgexport(h, 'report/figures/q1/q114-S0.eps');

h = figure('name','D');
imshow(mapD, [min(mapD(:)), max(mapD(:))])
hgexport(h, 'report/figures/q1/q114-D.eps');

h = figure('name','F');
imshow(mapF, [min(mapF(:)), max(mapF(:))])
hgexport(h, 'report/figures/q1/q114-F.eps');

h = figure('name','RESNORM');
imshow(mapRESNORM, [maskThreshMin, 300*maskThreshMin])
hgexport(h, 'report/figures/q1/q114-RESNORM.eps');


h = figure('name','fbDir');
quiver(fibDirXMap, fibDirYMap);
hgexport(h, 'report/figures/q1/q114-fbDir.eps');


end