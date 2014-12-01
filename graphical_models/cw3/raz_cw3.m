function [] = raz_cw3()

prob74()

end

function [] = prob74()

load('airplane.mat')

%DEMOMDP demo of solving Markov Decision Process on a grid
import brml.*
[Gx Gy] = size(U);  % two dimensional grid size
S = Gx*Gy; % number of states on grid
st = reshape(1:S,Gx,Gy); % assign each grid point a state

A = 5;  % number of action (decision) states
%[stay up down left right] = assign(1:A); % actions (decisions)
stay =1;
up =2;
down =3;
left =4;
right=5;
p = zeros(S,S,A); % initialise the transition p(xt|xtm,dtm) ie p(x(t)|x(t-1),d(t-1))

% make a deterministic transition matrix on a 2D grid:
for x = 1:Gx
	for y = 1:Gy
		p(st(x,y),st(x,y),stay)=1; % can stay in same state
		if validgridposition(x+1,y,Gx,Gy)
			p(st(x+1,y),st(x,y),right)=1;
		end
		if validgridposition(x-1,y,Gx,Gy)
			p(st(x-1,y),st(x,y),left)=1;
		end
		if validgridposition(x,y+1,Gx,Gy)
			p(st(x,y+1),st(x,y),up)=1;
		end
		if validgridposition(x,y-1,Gx,Gy)
			p(st(x,y-1),st(x,y),down)=1;
		end
	end
end
% define utilities

u = zeros(S,1);
for i=1:Gx
    for j=1:Gy
        u(st(i,j)) = U(i,j);
    end
end

gam = 0.95; % discount factor

%figure; imagesc(u); colorbar; title('utilities'); pause


[xt xtm dtm]=assign(1:3); % assign the variables x(t), x(t-1), d(t-1) to some numbers

% define the transition potentials p(x(t)|x(t-1),d(t-1))
tranpot=array([xt xtm dtm],p);
% setup the value potential v(x(t))
valpot=array(xt,ones(S,1)); % initial values

maxiterations=100; tol=0.001; % termination criteria

% Value Iteration:
oldvalue=valpot.table;
for valueloop=1:maxiterations
	valueloop
	tmppot = maxpot(sumpot(multpots([tranpot valpot]),xt),dtm);
	valpot.table = u + gam*tmppot.table; % Bellman's recursion
	if mean(abs(valpot.table-oldvalue))<tol; break; end % stop if converged
	oldvalue = valpot.table;
	%imagesc(reshape(valpot.table,Gx,Gy)); colorbar; drawnow
end

valpot.table

curr_pos = [1,13];
deltas = [[1,0], [0,1],[-1,0],[-1,-1],[0,0]]

best_path = []
while(curr_pos(1) ~= 8 && curr_pos(2) ~= 4)
    
    for delta = deltas
        new_pos = curr_pos + delta;
        if validgridposition(new_pos(1),new_pos(2),Gx,Gy)
            values = [values valpot.table(st(new_pos(1), new_pos(2)))];
        end
    end
    [~, max_index] = max(values);
    curr_pos = deltas(max_index) + curr_pos;
    best_path = [best_path curr_pos];
end

end