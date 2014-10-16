function [] = raz_probs()

prob27()
end

function [dec] = cliqueToDec(clique)

%binarray = zeros(10,1);
dec = 0;
for i=1:10
    if (ismember(i, clique))
        %binarray(i) = 1
        dec = dec + 2^(10-i);
    end
end

end

%% prob 2.6
function [] = prob26()
load('WikiAdjSmall.mat');

dist = graphallshortestpaths(A,'Directed',false);

hist = zeros(1000,1);

[width, height] = size(A);

for i=1:width
    for j=1:i
        if (dist(i,j) ~= inf && dist(i,j) ~= 0)
            hist(dist(i,j)) = hist(dist(i,j)) + 1;
        end
    end
end

hist(1:20)

end

%% prob 2.7

function [] = prob27()
C = load('cliques.mat');

C = C.cl;
len = length(C);

maxCliques = cell(100,1);
count = 1
for c1=1:len
    isContained = 0;
    for c2=1:len
        if(c1 ~= c2 && all(ismember(C{c1}, C{c2})))
            isContained = isContained + 1;
            break
        end
    end
    if(isContained == 0)
       maxCliques{count} = C{c1};
       count = count + 1;
    end
end

count = count - 1;

celldisp(maxCliques(1:count))

decCliques = zeros(count,1);
for i=1:count-1
    decCliques(i) = cliqueToDec(maxCliques{i});
end
    
decCliques



end


function [] = prob314()

[A B C]=assign(1:3)

variable(A).name='A'; variable(A).domain = {1,2};
variable(B).name='B'; variable(B).domain = {1,2,3};
variable(C).name='C'; variable(C).domain = {1,2};

pot(A) = array;
pot(A).variables=[1,2];
pot{A}.table(1) = 1/2;
pot{A}.table(2) = 1/2;

pot(B).variables=[1 2 3]; % define array below using this variable order
tmptable(yes, yes)=0.05; tmptable(yes, no)=0.01;  
tmptable(no,:)=1-tmptable(yes,:); % due to normalisation
pot(tub).table=tmptable; 

end