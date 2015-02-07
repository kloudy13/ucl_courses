function q114_generate_maps()

load('q114.mat')

maskThreshMin = 1*10^8;
maskThreshMax = 4*10^9;
mask = (mapRESNORM > maskThreshMin) .* (mapRESNORM < maskThreshMax);

%figure('name','Fbef');
%imshow(mapF, [min(mapF(:)), max(mapF(:))])

%mapRESNORM = mapRESNORM .* mask;
%mapD = mapD .* mask;
mapF = mapF .* mask;
%mapS0 = mapS0 .* mask;

%figure('name','Fafter');
%imshow(mapF, [min(mapF(:)), max(mapF(:))])

fibDirXMap = mapF .* cos(mapPhi) .* sin(mapTheta);
fibDirYMap = mapF .* sin(mapPhi) .* sin(mapTheta);
fibDirXMap = fibDirXMap .* mask;
fibDirYMap = fibDirYMap .* mask;

h = figure('name','S0');
imshow(mapS0, [min(mapS0(:)), max(mapS0(:))])
%hgexport(h, 'report/figures/q1/q114-S0.eps');
saveTightFigure(h,'report/figures/q1/q114-S0.eps')

h2 = figure('name','D');
imshow(mapD, [0.5e-03, 5e-03])
saveTightFigure(h2, 'report/figures/q1/q114-D.eps');

h3 = figure('name','F');
imshow(mapF, [min(mapF(:)), max(mapF(:))])
saveTightFigure(h3, 'report/figures/q1/q114-F.eps');

h4 = figure('name','RESNORM');
imshow(mapRESNORM, [maskThreshMin, 100*maskThreshMin])
saveTightFigure(h4, 'report/figures/q1/q114-RESNORM.eps');


h5 = figure('name','fbDir');
quiver(fibDirXMap, fibDirYMap);
saveTightFigure(h5, 'report/figures/q1/q114-fbDir.eps');


end