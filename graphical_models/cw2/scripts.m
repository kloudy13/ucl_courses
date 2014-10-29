function scripts()




prob120()

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