function scripts()


%prob57()


prob67()
%prob120()
%pro69() - checked, gives correct result

%pro69cv2()

end


function prob120()

len = 10;
grid = zeros(len, len); % in each cell (i,j) contains how many possible scenarios exist where it is occupied by any of the boats

boatLen = 5;

illegal_locations = [1, 10; 2, 2;3, 8; 4, 4; 5, 6; 6, 5; 7, 4; 7, 7; 9, 2; 9, 9]';

% a few tests
testPoints = [ones(boatLen, 1) * 1, ones(boatLen, 1) * 6 + (0:boatLen-1)']'; 
assert(boat_illegal_location(testPoints, illegal_locations ) == 1);

testPoints = [ones(boatLen, 1) * 1, ones(boatLen, 1) * 5 + (0:boatLen-1)']'; 
assert(boat_illegal_location(testPoints, illegal_locations ) == 0);


% horizontal boat at position (iH,jH) occupies cells (iH, jH), (iH, jH+1), (iH,
% jH+2), (iH, jH+3) ...
% vertical boat at position (iV,jV) occupies cells (iV,jV), (iV+1,jV),
% (iV+2,jV) ...

% fix location of the horizontal ship
for iH=1:len
    for jH=1:(len-boatLen+1)
        % find cells the boat is occupying
        hBoatPoints = [ones(boatLen, 1) * iH, ones(boatLen, 1) * jH + (0:boatLen-1)']';  
        
        % make sure boat is not on illegal location
        if (boat_illegal_location(hBoatPoints, illegal_locations))
           continue; 
        end    
        
        % fix location of the vertical ship
        for iV=1:(len-boatLen+1)
            for jV=1:len
                % find cells the boat is occupying
                vBoatPoints = [ones(boatLen, 1) * iV  + (0:boatLen-1)', ones(boatLen, 1) * jV]'; 
                        
                % make sure boat is not on illegal location
                if (boat_illegal_location(vBoatPoints, illegal_locations))
                   continue; 
                end    
                
                % make sure boats don't collide
                if (boat_illegal_location(vBoatPoints, hBoatPoints))
                    %vBoatPoints
                    %hBoatPoints
                   continue; 
                end

                % update grid counters
                %[hBoatPoints, vBoatPoints]
                grid = update_grid(grid, [vBoatPoints, hBoatPoints]);
            end
        end
    end
end


grid = grid ./ sum(grid(:));

% multiply everything by 10 because there are 10 pixels that are actually
% occupied on the board (5 pixels for each ship)
grid = grid .* 10;
grid

max = 0;
maxI = 1; 
maxJ = 1;

for i=1:len
    for j=1:len
        if(grid(i,j) > max)
           maxI = i;
           maxJ = j;
           max = grid(i,j); 
        end
    end
end

max
[maxI, maxJ]

end

function is_on_illegal_loc = boat_illegal_location(boatPoints, illegal_locations)

is_on_illegal_loc = 0;

for point=boatPoints
    for loc=illegal_locations
        if (any(point - loc) == 0)
            is_on_illegal_loc = 1;
            return;
        end
    end
    
end

end

function grid = update_grid(grid, boatPoints)

for p=boatPoints
    grid(p(1),p(2)) = grid(p(1),p(2)) + 1;
end

end

function [] = prob57()

format long

load('banana.mat');

T = length(x);

%DEMOSUMPROD Sum-Product algorithm test :
import brml.*
% Variable order is arbitary

xVars=1:T;
yVars=T+1:2*T;
hVars=2*T+1:3*T;

xstates=1:4;
ystates=1:4;
hstates=1:5;

xGh = 1; 
htGhtm = 2;
htm = 3;

% maxHt(t, hState, :) stores [max, argmax(ht)] at iteration t 
maxHt = zeros(T, hstates, 2); % the third dimension is two because we store both 

% initialise the potential tables p(x|h), p(ht|ht-1), p(h1) 
pot{xGh}.table=pxgh;
pot{htGhtm}.table=phtghtm;
pot{htm}.table=ph1;



for t=2:T
    xtState = charToInt(x(t-1));
    
    pot{xGh}.variables=[xVars(t-1) hVars(t-1)];
    pot{htGhtm}.variables=[hVars(t) hVars(t-1)];
    pot{htm}.variables=[hVars(t-1)];

    pot=setpotclass(pot,'array');
    
    % multiply the potentials jointpot = p(x|h) * p(ht|ht-1) * p(ht)
    jointpot = multpots(pot);

    %potHtHtmGx = condpot(jointpot, [hVars(t) hVars(t-1)], xVars(t-1) )

    %potHtHtm = setpot(potHtHtmGx ,xVars(t-1), xtState)
    % set xt to the observed value and get p(ht, ht-1, x = observed_x)
    potHtHtm = setpot(jointpot, xVars(t-1), xtState);

    %potHtgHtm.table

    % for each state of ht calculate the argmax_{ht-1} [ p(ht, ht-1) ]
    for htFixed=hstates
        % fix state of ht to htState
        htArray = setpot(potHtHtm, hVars(t), htFixed).table;
        [maxHt(t-1, htFixed, 1), maxHt(t-1, htFixed, 2)] = max(htArray);
    end

    % make the max a prior on ht and normalise it
    pHt = maxHt(t-1, :, 1)
    pHt = pHt ./ sum(pHt)
    
    pot{htm}.table=pHt;

end

maxHt

t = T+1;

xtState = charToInt(x(t-1));

pot2{xGh}.variables=[xVars(t-1) hVars(t-1)];
pot2{htm}.variables=[hVars(t-1)];
pot2{xGh}.table=pxgh;
pot2{htm}.table=pHt;

pot2=setpotclass(pot2,'array');

% multiply the potentials jointpot = p(x|h) * p(ht|ht-1) * p(ht)
jointpot = multpots(pot2);

%potHtHtmGx = condpot(jointpot, [hVars(t) hVars(t-1)], xVars(t-1) )

%potHtHtm = setpot(potHtHtmGx ,xVars(t-1), xtState)
% set xt to the observed value and get p(ht, ht-1, x = observed_x)
potHT = setpot(jointpot, xVars(t-1), xtState);

hStar = zeros(T, 1);

[~, hStar(T)] = max(potHT.table);

%% now propagate back all the max values
for i=fliplr(1:T-1)
    maxHt(i, hStar(i+1), 2)
    hStar(i) = maxHt(i, hStar(i+1), 2);
end

hStar

% calc argmax p(y|h*) = [p(y_1|h_1*), p(y_2|h_2*) ... ] proof given in the report 
yStar = zeros(T,1)
for t=1:T
    [~, yStar(t)] = max(pygh(:, hStar(t)));
end

yStarStr = [];

letterMap = 'ACGT';

for i=1:T
   yStarStr(i) = letterMap(yStar(i)); 
end

yStarStr = char(yStarStr)

konstantinYStar = 'CTTGACTTGACTGACTGACTGACTGACCTGATTTTTTGACTGAGACTGACTGTTTTTTTTCTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA';

assert(strcmp(yStarStr, konstantinYStar))
end


function i = charToInt(ch)
i = 0;
if ch == 'A'
    i = 1;
end
if ch == 'C'
    i = 2;
end
if ch == 'G'
    i = 3;
end
if ch == 'T'
    i = 4;
end

if i == 0
   throw(Exception);    
end
    
end



function prob67()
    n = 10;

    msg = ones(2^n,1);
    psyMatrix = buildPsyMatrix(n);
    psyMatrixFinal = buildPsyMatrixFinal(n);
        
    save('69_psy_matrices.mat', 'psyMatrix', 'psyMatrixFinal')
    %load('69_psy_matrices.mat')
    
    for iter = 2:9;
        msg(1:20)
        msg = psyMatrix * msg;
    end
    
    msg = psyMatrixFinal * msg;
    
    Z = sum(msg(:))
    logZ = log(Z)
    
end

function msgImtI = prob67aux(msgImtIm, iter, n)

msgImtI = zeros(2^n,1);
psy= zeros(2^n,1);

for xI=0:2^n-1
    bigSum = 0;
    vI = binary2vector(xI, n); % binary values of the small xIs as a vector of n elements
    
    for xIm=0:2^n-1
        vIm = binary2vector(xIm, n);
        if (iter<n)
            psy(xIm+1) = exp(sum((vI - vIm) .^2) + sum((vIm(1:n-1) - vIm(2:end)) .^ 2));
        else
            psy(xIm+1) = exp(sum((vIm(1:n-1) - vIm(2:end)) .^ 2));
        end
    end
    
    msgImtI(xI+1) = sum(psy .* msgImtIm);
end


end

function psyMatrix = buildPsyMatrix(n)

psyMatrix = zeros(2^n,2^n);

for xI=0:2^n-1
    vI = binary2vector(xI, n); % binary values of the small xIs as a vector of n elements
    xI
    for xIm=0:2^n-1
        vIm = binary2vector(xIm, n);

        psyMatrix(xI+1, xIm+1) = exp(sum((vI - vIm) .^2) + sum((vIm(1:n-1) - vIm(2:end)) .^ 2));

    end
   
end
end

function psyMatrix = buildPsyMatrixFinal(n)

psyMatrix = zeros(2^n,2^n);

for xI=0:2^n-1
    vI = binary2vector(xI, n); % binary values of the small xIs as a vector of n elements
    xI
    for xIm=0:2^n-1
        vIm = binary2vector(xIm, n);

        psyMatrix(xI+1, xIm+1) = exp(sum((vI - vIm) .^2) + sum((vIm(1:n-1) - vIm(2:end)) .^ 2) + sum((vI(1:n-1) - vI(2:end)) .^ 2));

    end
   
end
end

function [] = pro69()

format long
import brml.*
load('diseaseNet')

% code snippets taken from demoJTree.m

pot=str2cell(setpotclass(pot,'array')); % convert to cell array 
%pot=setpotclass(pot,'array'); % convert to cell array 


[jtpot jtsep infostruct]=jtree(pot); % setup the Junction Tree

[jtpot jtsep logZ]=absorption(jtpot,jtsep,infostruct); % do full round of absorption

nrSeps = 40;
pS = zeros(nrSeps, 1);

for var=1:nrSeps
    jtpotnum = whichpot(jtpot,var+20,1); % find a single JT potential that contains dys
    tmpTable=table(sumpot(jtpot(jtpotnum),var+20,0)); % sum over everything but dys
    pS(var) = double(tmpTable(1));
end

pS


%    0.441834190079851
%    0.456675263152281
%    0.441405071163811
%    0.491274625876646
%    0.493892068515409
%    0.657483115182047
%    0.504562541586642
%    0.268693391658929
%    0.649077145409492
%    0.490740180661327
%    0.422551897884589
%    0.429096175062460
%    0.545020871745821
%    0.632958925914310
%    0.429540436814508
%    0.458794090333332
%    0.427558815937816
%    0.404255167991632
%    0.582092524680478
%    0.589590513057441
%    0.761270142889589
%    0.695587782104791
%    0.508701848643658
%    0.419961576354167
%    0.351942467409690
%    0.389611106637832
%    0.325972664362974
%    0.469624102888099
%    0.522867791671208
%    0.717311545958481
%    0.524199126536054
%    0.353703968717321
%    0.512679203456453
%    0.529404158281906
%    0.385750990816217
%    0.489094812970005
%    0.633595387016685
%    0.589604428141779
%    0.423164229556478
%    0.528234316630121

pS2 = zeros(nrSeps, 1);

for var=1:nrSeps
    potVar = pot{var+20};
    %normpot = potVar ./ sum(potVar(:));
    %normpot = normp(potVar)
    parents = pot{var+20}.variables(2:end);
    jointPot{1}=potVar;
    jointPot{2}=pot{parents(1)};
    jointPot{3}=pot{parents(2)};
    jointPot{4}=pot{parents(3)};
    potVar = multpots(jointPot);
%     if(var == 1)
%        potVar 
%     end
    tmpTable=table(sumpot(potVar,var+20,0)); % sum over everything but dys
    pS2(var) = double(tmpTable(1));
end

pS2

for i=1:length(pS)
    assert(abs(pS(i) - pS2(i)) < 0.000001)
end

end

function [] = pro69c()

format long
import brml.*
load('diseaseNet')

% code snippets taken from demoJTree.m

symptoms = [1,1,1,1,1,2,2,2,2,2];

pot=str2cell(setpotclass(pot,'array')); % convert to cell array 
%pot=setpotclass(pot,'array'); % convert to cell array 

nrDiseases = 20;
pDigS = zeros(nrDiseases, 2);

for i=1:nrDiseases
   %display(strcat('Iteration ', str(i)))
   i
   for di=0:1
       sum = 0;
       for d=0:2^(nrDiseases-1) - 1
          module_size = 100; 
          if (mod(d, module_size) == 0)
             fprintf('progress %d\n', d / module_size);
          end
          binD = binary2vector(d, nrDiseases-1);
          binD = [binD(1:i-1), di, binD(i:end)];
          assert(length(binD) == nrDiseases);
          prod = 1;
          for sIndex=1:length(symptoms)
              parentDiseases = pot{sIndex + 20}.variables(2:end);
              prod = prod * pot{sIndex + 20}.table(symptoms(sIndex), binD(parentDiseases(1))+1, binD(parentDiseases(2))+1, binD(parentDiseases(3))+1);
          end
          for dIndex=1:nrDiseases
              prod = prod * pot{dIndex}.table(binD(dIndex)+1);
          end
          sum = sum + prod;

       end
       pDigS(i, di+1) = sum;
   end
    
end

pDigS

end


function out = binary2vector(data,nBits)

powOf2 = 2.^[0:nBits-1];

%# do a tiny bit of error-checking
if data > sum(powOf2)
   error('not enough bits to represent the data')
end

out = false(1,nBits);

ct = nBits;

while data>0
if data >= powOf2(ct)
data = data-powOf2(ct);
out(ct) = true;
end
ct = ct - 1;
end

end


function [] = pro69cv2()

clear all; clc;
format long
import brml.*
load('diseaseNet')

% code snippets taken from demoJTree.m

pot=str2cell(setpotclass(pot,'array')); % convert to cell array 
%pot=setpotclass(pot,'array'); % convert to cell array 

%pot = pot(1:30);

symptoms = [1,1,1,1,1,2,2,2,2,2];

displacement = 20;


for s=1:10
   pot{s+displacement}.table
   display('-------------------------')
   pot{s+displacement} = setpot(pot{s+displacement}, s+displacement, symptoms(s));
   pot{s+displacement}.table
   sum(pot{s+displacement}.table(:))
end

%    0.029775615005615   0.970224384994385
%    0.000225780279841   0.000365638819622
%    0.000564352736370   0.000027066363093
%    0.000234583097360   0.000356836002103
%    0.000293620346471   0.000297798752993
%    0.000257358619973   0.000334060479490
%    0.000110883531100   0.000480535568364
%    0.000414692910648   0.000176726188816
%    0.000025505923381   0.000565913176082
%    0.000360950480353   0.000230468619111
%    0.000169927892654   0.000421491206809
%    0.000289696707173   0.000301722392291
%    0.000532040521213   0.000059378578250
%    0.000366422335790   0.000224996763674
%    0.920475727578874   0.079524272421126
%    0.000417598757918   0.000173820341546
%    0.000119021441593   0.000472397657870
%    0.000537300889925   0.000054118209539
%    0.000511558152920   0.000079860946543
%    0.000522772520397   0.000068646579067

%    0.029775615005615   0.970224384994385
%    0.381760210392464   0.618239789607536
%    0.954234885011638   0.045765114988362
%    0.396644439743479   0.603355560256521
%    0.496467474143329   0.503532525856671
%    0.435154394246893   0.564845605753107
%    0.187487234011107   0.812512765988893
%    0.701182817774887   0.298817182225113
%    0.043126648098117   0.956873351901883
%    0.610312518956490   0.389687481043509
%    0.287322294474943   0.712677705525057
%    0.489833195166369   0.510166804833631
%    0.899599829792365   0.100400170207635
%    0.619564596615449   0.380435403384551
%    0.920475727578874   0.079524272421126
%    0.706096164794739   0.293903835205261
%    0.201247206424932   0.798752793575068
%    0.908494315473291   0.091505684526709
%    0.864967251453648   0.135032748546352
%    0.883929046037277   0.116070953962723


[jtpot jtsep infostruct]=jtree(pot); % setup the Junction Tree

[jtpot jtsep logZ]=absorption(jtpot,jtsep,infostruct); % do full round of absorption

nrDiseases = 20;
pDgS = zeros(nrDiseases, 2);

for var=1:nrDiseases
    jtpotnum = whichpot(jtpot,var,1); % find a single JT potential that contains dys
    tmpTable=table(sumpot(jtpot(jtpotnum),var,0)); % sum over everything but dys
    pDgS(var,:) = double(tmpTable);
end

% wrong, need to take the values of the symptoms into account!!!

pDgS

for i=1:length(pDgS)
    pDgS(i,:)  = pDgS(i,:) ./ sum(pDgS(i,:));
end
pDgS

end

function [] = prob513()

load('SimoHurrta')

% p(S_{T} | S_{T-1})
pStgStm = 

end
