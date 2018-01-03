% cellular model that will ultimately be a spatially distributed albedo
% scheme

% 1. set up random grid
% 2. update by growing blooms probabalistically
% 3. use values in grid to call appropriate snicar
% 4. iteratively update second grid with albedo values

% %% TODOS
    % Deal with edge effects (-1s and +1s can push the indexing out of grid
    % range, currently counter starts at 2 and ends at i,j -1, but this
    % leads to a border that doesn't update)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on

% %% Set up grid
meanlist = []; % list to append mean grid value to for surface-class grid

gridsize = 10000; % total grid area (no. cells)
gridx = 100; % length of x-axis
gridy= 100; % length of y-axis
alg_frac = 1; % percentage of algal coverage at start of experiment (all initialise as light algae: class 1)
non_alg_frac = gridsize-alg_frac; % residual = non-algal, assumed clean

% create grid by randomly spacing 1's and 2's in a 10x10 grid
alg_pixels = gridsize*(alg_frac/100) 
non_alg_pixels = gridsize-alg_pixels

grid = reshape([zeros(non_alg_pixels,1) ; ones(alg_pixels,1)],gridx,gridy) ;
grid(:) = grid(randperm(numel(grid))) ;

for counter = 1:1:12
    rand = randi(100,1); % return one pseudorandom integer between 1-100
    for i = 2:1:gridx-1  %%%% NEED TO DEAL WITH EDGES BETTER!!!!
        for j = 2:1:gridy-1 %%%% NEED TO DEAL WITH EDGES BETTER!!!!        
            
            if grid(i,j) == 0   % if cell = 0, 10% chance of bloom initiation, otherwise stay clean
                if rand == 1
                    grid(i,j) = grid(i,j)+1;
                else
                    grid(i,j) = grid(i,j);
                end
                
            elseif grid(i,j) < 10 % if cell is between 1 and 10, 40% chance of increasing in situ
                if rand < 30
                    grid(i,j) = grid(i,j)+1;
                elseif rand > 30 && rand < 80 % 40 % chance of bloom spreading, not darkening
                    rand1 = randi(8,1); % each neighbour has equal chance of being colonised
                    if rand1 == 1
                        grid(i,j+1) = grid(i,j+1)+1;
                    elseif rand1 == 2
                        grid(i,j-1) = grid(i,j-1)+1;
                    elseif rand1 == 3
                        grid(i+1,j+1) = grid(i+1,j+1)+1;
                    elseif rand1 == 4
                        grid(i+1,j-1) = grid(i+1,j-1)+1;
                    elseif rand1 == 5
                        grid(i+1,j) = grid(i+1,j)+1;
                    elseif rand1 == 6
                        grid(i-1,j) = grid(i-1,j)+1;
                    elseif rand1 == 7
                        grid(i-1,j-1) = grid(i-1,j-1)+1
                    elseif rand1 == 8
                        grid(i-1,j+1) = grid(i,j+1)+1
                    end

            elseif grid(i,j) == 10 % if cell value is 10 bloom can only grow by spreading
                if i < gridx && j < gridy && i >1 && j >1
                    if rand1 == 1
                        grid(i,j+1) = grid(i,j+1)+1;
                    elseif rand1 == 2
                        grid(i,j-1) = grid(i,j-1)+1;
                    elseif rand1 == 3
                        grid(i+1,j+1) = grid(i+1,j+1)+1;
                    elseif rand1 == 4
                        grid(i+1,j-1) = grid(i+1,j-1)+1;
                    elseif rand1 == 5
                        grid(i+1,j) = grid(i+1,j)+1;
                    elseif rand1 == 6
                        grid(i-1,j) = grid(i-1,j)+1;
                    elseif rand1 == 7
                        grid(i-1,j-1) = grid(i-1,j-1)+1;
                    elseif rand1 == 8
                        grid(i-1,j+1) = grid(i,j+1)+1;
                    end
                else
                    grid(i,j) = grid(i,j);
                end
                end
            end
        end
    end
                
    meanlist(counter) = mean2(grid);
    surf(grid);
    view(2);
    cbar = colorbar;
    cbar.Limits=[0,10];
    caxis([0,10]);
    drawnow
    pause(0.1)
    

    
    % calls to SNICAR should go in here: populate second grid with albedos
    
    
end
