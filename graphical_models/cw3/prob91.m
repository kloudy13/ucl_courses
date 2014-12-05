function prob91()

import brml.*

load('printer.mat')

[nr_nodes, nr_visits] = size(x);
x = x-ones(nr_nodes, nr_visits);

[fuse drum toner paper roller burn quality wrinkled mult jam] = assign(1:nr_nodes);

pots = cell(nr_nodes,1);


deps = cell(nr_nodes,1);
deps{fuse} = [];
deps{drum} = [];
deps{toner} = [];
deps{paper} = [];
deps{roller} = [];
deps{burn} = [fuse];
deps{quality} = [drum, toner, paper];
deps{wrinkled} = [fuse, paper];
deps{mult} = [paper, roller];
deps{jam} = [fuse, roller];

for node=1:nr_nodes
    % filter x 
    parents = deps{node};
    nr_parents = length(parents);
    new_x = x;
    if nr_parents > 0
        potTable = zeros(2 * ones(1,nr_parents+1));
        for state=0:2^nr_parents - 1
            binState = binary2vector(state,nr_parents);
            indicesToKeep = ones(1, nr_visits);
            for pI=1:nr_parents
                parent = parents(pI);
                indicesToKeep = indicesToKeep .* (x(parent,:) == binState(pI));
            end 
            new_x = x(:,indicesToKeep);
            
            if(nr_parents == 1)
                potTable(:,binState(1)) = new_x;
            end
            if(nr_parents == 2)
                potTable(:,binState(1),binState(2)) = new_x;
            end
            if(nr_parents == 3)
                potTable(:,binState(1),binState(2),binState(3)) = new_x;
            end
        end
        pots{node} = array([node parents], getPot(potTable, node)); 
    end
    
    if nr_parents == 0
        pots{node} = array([node], getPot(x, node));
    end
end



end

function pot = getPot(x, node)
[nr_nodes, nr_visits] = size(x);
pot = [sum(x(node,:)) nr_visits - sum(x(node,:))];
pot = pot ./ sum(pot);
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
