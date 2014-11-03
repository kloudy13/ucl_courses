function scripts()


%prob57()


%prob67()
%prob120()
pro69()

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


msgImtI = ones(2^n,1);



i = 2;


end

function msgImtI = prob67aux(msgImtI, iter, n)


for xI=1:2^n
    bigSum = 0;
    vI = de2bi(xI, n); % binary values of the small xIs as a vector of n elements
    for xIm=1:2^n
        vIm = de2bi(xIm, n); 
        
        % go vertical
        prod = 1;
        for xVert=1:n-1
           %prod = prod * 
        end
        
        % go horizontal
    end
end

end

function [] = pro69()

import brml.*
load('diseaseNet')

% code snippets taken from demoJTree.m

pot=str2cell(setpotclass(pot,'array')); % convert to cell array 

[jtpot jtsep infostruct]=jtree(pot); % setup the Junction Tree

[jtpot jtsep logZ]=absorption(jtpot,jtsep,infostruct); % do full round of absorption

%figure; drawNet(dag(pot),variable); title('Belief Net');
%figure; drawJTree(infostruct,variable); title('Junction Tree (separators not shown)');

% for pot=jtpot
%     for var=pot.variables
%        if pS(var) == 0
%           pS(var) =  
%        end
%     end
% end


nrSeps = 40;
pS = zeros(nrSeps, 1);

for var=1:nrSeps
    jtpotnum = whichpot(jtpot,var+20,1); % find a single JT potential that contains dys
    tmpTable=table(sumpot(jtpot(jtpotnum),var+20,0)); % sum over everything but dys
    pS(var) = double(tmpTable(1));
end

pS

pS2 = zeros(nrSeps, 1);



for var=1:nrSeps
    potVar = pot{var+20};
    %normpot = potVar ./ sum(potVar(:));
    %normpot = normp(potVar)
    tmpTable=table(sumpot(potVar,var+20,0)); % sum over everything but dys
    pS2(var) = double(tmpTable(1));
end

pS2 = pS2 ./ 8




end