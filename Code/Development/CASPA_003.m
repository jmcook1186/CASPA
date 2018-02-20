% 
% Radiative Transfer model distributed over space and time

% %% OVERVIEW
% This model simulates a snow surface with physical/optical properties
% defined by the user. This surface is represented as a cellular automaton,
% i.e. an x,y grid composed of discrete cells that update based upon certain
% functions as the model advances through time.
% There are various modes of operation. The simplest is to simulate an
% algal bloom that grows from an initial configuration defined by the user.
% The user defines the % coverage at t=0. At each timestep, algal growth is
% simulated using a probabalistic function, where growth in situ has a
% given likelihood, spread of the algae into neighbouring cells has a given
% likelihood, and spontaneous initiation of a new bloom has a give
% likelihood. Growth is represented by increasing the cell value - i.e. 0
% represents clean ice, 1-10 represent increasing mixing ratio of algae in
% the upper surface. A carrying capacity is simulated wherein the cell
% value cannot exceed a value of 10. Once the cell value reaches 10 the
% bloom can only grow by spreading. If all the neighbouring cells are also
% at carrying capacity, the algae will not grow.

% At each timestep, the grid is thereby updated. The value of each cell is
% associated with a call to a specific instance of the radiative transfer
% scheme BioSNICAR. These instances are pre-coded as driver routines
% saved in the working directory. In this default version, each instance is
% identical except for the mixing ratio of algae in the upper 3mm of ice.
% The broadband and spectral albedo of each pixel in the grid is thereby 
% modelled using BioSNICAR and a spatially-integrated albedo computed 
% using a grid-mean.

% Also, in the default version the algal growth is uninterupted, so left
% long enough the algae will have even coverage at carrying capacity.
% However, there is the option to interrupt the growth, simulating a
% rainfall event for example. Could potentially simulate nutrient
% limitation by increasin or decreasing the doubling times for the blooms.
% A nutrient nutrient pulse could thereby be simulated.

% Possibility to model indirect effects: could increase grain size along
% with algal biomass, according to field observations (measured biomass vs
% depth of 1030nm absorption feature). 

% %% TODOS
    
    % IDEA: make each instance of snicar_alg precisely double the biomass
    % of the previous (i.e. 2 = 2x1, 3 = 2x2, 4 = 2x3...). This way we can
    % relate time to real algal doubling times - max 3 days, min 7days.
    % Then, each timestep = 3-7 days. Season = 3 months 
    % = 12 timesteps for 7 day doubling time. For 3 days doubling time a 3
    % month season is simulated as 28 timesteps (28 x 3 days = 3 months)
    
    % 1. Deal with edge effects (-1s and +1s can push the indexing out of grid
    % range, currently counter starts at 2 and ends at i,j -1, but this
    % leads to a border that doesn't update)
    
    % 2. plot(meangrid) after runs and analyse temporal evolution
    
    % 3. Could add an additional randi() dice roll to randomly reduce value
    % of cells, simulating rainfall washaway. Alternatively, a sudden dust
    % deposition event...? These don't have to be random, they could start
    % on a particular ticker value (e.g. at 20th timestep: rain).
        
    % 4. Can the initial grid be set using UAV por satellite imagery?
    
    % 5. How should I deal with dust? constant background? 
    
    % 6. Possibility of calculating radiative forcing in each SNICAR run,
    % then calculating energy abalnce and volume of melt,. and linking that
    % to change in grain size? Iterate through and incorporate indirect
    % albedo evolution.

    
%%%  VERSION 3 EDITS:
% 1) No longer runs SNICAR per cell. Instead, all the instances of SNICAR used
% in the model are run once at the beginning and the spectral and broadband
% data appended to lists. These are then used as a lookup table for the
% cellular model. Of course, this has drastically reduced the run time.
%
% 2) Bug fix in final loop - previously able to increase cell values above
% 10. This was because when an individual cell reached carrying capacity
% the algae grew by spreading but was able to spread into cells whose
% values > 10. This is now fixed so no cell can exceed the carrying
% capacity.
%
% 3) Condense loops - now single loop for if rand>50 OR cell value 
% =10 (previously two sequential loops)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a figure for plotting the final albedo grid
figure;

% set up empty lists

BBAlist = []; % list to append mean grid value to for albedo grid

%%% USER DEFINED VARIABLES: SET HERE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up time ticker (time_tot = total no of timesteps in run)
time_tot = 200; 
timestep = 1;

% %% Set up grid
gridsize = 10000; % total grid area (no. cells)
gridx = 100; % length of x-axis
gridy= 100; % length of y-axis
alg_frac = 1; % percentage of algal coverage at start of experiment (all initialise as light algae: class 1)
non_alg_frac = gridsize-alg_frac; % residual = non-algal, assumed clean

% turn % coverage into actual number of pixels
alg_pixels = gridsize*(alg_frac/100); 
non_alg_pixels = gridsize-alg_pixels;

% create grid and randomly place 0's and 1's according to user defined
% coverage
grid = reshape([zeros(non_alg_pixels,1) ; ones(alg_pixels,1)],gridx,gridy) ;
grid(:) = grid(randperm(numel(grid))) ;

% start stepping through time
for counter = 1:timestep:time_tot    
    
    rand = randi(100,1); % return one pseudorandom integer between 1-100 (each integer has 1% chance of selection)
    
    for i = 2:1:gridx-1  %%%% NEED TO DEAL WITH EDGES BETTER!!!!
        for j = 2:1:gridy-1 %%%% NEED TO DEAL WITH EDGES BETTER!!!!            
      
            if grid(i,j) == 0   % if cell = 0, stay clean
                grid(i,j) = grid(i,j);
                       
            elseif grid(i,j) > 0 && grid(i,j) < 10 % if cell is between 1 and 9, 50% chance of increasing in situ
                if rand < 50
                    grid(i,j) = grid(i,j)+1;
                    
                elseif rand > 50 || grid(i,j) == 10 % if rand > 50 OR cell value = 10, bloom spreads, not darkening in situ
                    rand1 = randi(8,1); % each neighbour has equal chance of being colonised
            
                    if i < gridx && j < gridy && i >1 && j >1
                        if rand1 == 1
                            if grid(i,j+1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i,j+1)= grid(i,j+1);
                            else
                                grid(i,j+1) = grid(i,j+1)+1;
                            end
                        elseif rand1 == 2
                            if grid(i,j-1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i,j-1)= grid(i,j-1);
                            else
                                grid(i,j-1) = grid(i,j-1)+1;
                            end
                        elseif rand1 == 3
                            if grid(i+1,j+1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i+1,j+1)= grid(i+1,j+1);
                            else
                                grid(i+1,j+1) = grid(i,j+1)+1;
                            end
                        elseif rand1 == 4
                            if grid(i+1,j-1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i+1,j-1)= grid(i+1,j-1);
                            else
                                grid(i+1,j-1) = grid(i+1,j-1)+1;
                            end
                        elseif rand1 == 5
                            if grid(i+1,j) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i+1,j)= grid(i+1,j);
                            else
                                grid(i+1,j) = grid(i+1,j)+1;
                            end
                        elseif rand1 == 6
                            if grid(i-1,j) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i-1,j)= grid(i-1,j);
                            else
                                grid(i-1,j) = grid(i-1,j)+1;
                            end
                        elseif rand1 == 7
                            if grid(i-1,j-1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i-1,j-1)= grid(i-1,j-1);
                            else
                                grid(i-1,j-1) = grid(i-1,j-1)+1;
                            end
                        elseif rand1 == 8
                            if grid(i-1,j+1) >=10 % check if the random cell is already at carrying capacity, if it is, no change
                                grid(i-1,j+1)= grid(i-1,j+1);
                            else
                                grid(i-1,j+1) = grid(i-1,j+1)+1;
                            end
                        end
                    else
                        grid(i,j) = grid(i,j);
                    end
                end
            end
        end
    end

% %% Plotting for initial 'surface class' grid disabled. Uncomment to show plot.

%     meanlist(counter) = mean2(grid);
%     surf(grid);
%     view(2);
%     cbar = colorbar;
%     cbar.Limits=[0,10];
%     caxis([0,10]);
%     drawnow
%     pause(0.2)

%%%%%%% CUT MODEL HERE IF SURFACE CLASS GRID ONLY REQUIRED %%%%%%%%%%%%%%%
                     %%%%%%%%%%%%%%%%%%%%%%
                     
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% RUN SNICAR AND UPDATE NEW GRID WITH ALBEDO VALUES %%%%%%%%%%%
    
grid2 = grid; % copy surface class grid for populating with albedos

%%% RUN ALL POSSIBLE PERMUTATIONS OF SNICAR BEFORE LOOPING THROUGH GRID %%%%%
[snicar_clean_BBA,snicar_clean_spectral]=SNICAR_driver_clean_CA;
[snicar_alg1_BBA,snicar_alg1_spectral]=SNICAR_driver_algae_CA_1;
[snicar_alg2_BBA,snicar_alg2_spectral]=SNICAR_driver_algae_CA_2;
[snicar_alg3_BBA,snicar_alg3_spectral]=SNICAR_driver_algae_CA_3;
[snicar_alg4_BBA,snicar_alg4_spectral]=SNICAR_driver_algae_CA_4;
[snicar_alg5_BBA,snicar_alg5_spectral]=SNICAR_driver_algae_CA_5;
[snicar_alg6_BBA,snicar_alg6_spectral]=SNICAR_driver_algae_CA_6;
[snicar_alg7_BBA,snicar_alg7_spectral]=SNICAR_driver_algae_CA_7;
[snicar_alg8_BBA,snicar_alg8_spectral]=SNICAR_driver_algae_CA_8;
[snicar_alg9_BBA,snicar_alg9_spectral]=SNICAR_driver_algae_CA_9;
[snicar_alg10_BBA,snicar_alg10_spectral]=SNICAR_driver_algae_CA_10;
        
%%% Loop through grid and replace surface class with albedo from SNICAR lookup
% library populated in loop above.
        
    for i = 1:1:gridx  % Loop through all cells
        for j = 1:1:gridy
            if grid(i,j) == 0  % if class = 0, run clean ice snicar, return BBA and spectral albedo 
                grid2(i,j) = snicar_clean_BBA; % add BBA to appropriate cell in grid2
            elseif grid2(i,j) ==1 % if cell is not 0, run appropriate 'algal' version of snicar
                grid2(i,j) = snicar_alg1_BBA;
            elseif grid(i,j) ==2
                grid2(i,j) = snicar_alg2_BBA;
            elseif grid(i,j) ==3
                grid2(i,j) = snicar_alg3_BBA;
            elseif grid(i,j) ==4
                grid2(i,j) = snicar_alg4_BBA;
            elseif grid(i,j) ==5
                grid2(i,j) = snicar_alg5_BBA;
            elseif grid(i,j) ==6
                grid2(i,j) = snicar_alg6_BBA;
            elseif grid(i,j) ==7
                grid2(i,j) = snicar_alg7_BBA;                
            elseif grid(i,j) ==8
                grid2(i,j) = snicar_alg8_BBA;
            elseif grid(i,j) ==9
                grid2(i,j) = snicar_alg9_BBA;
            elseif grid(i,j) ==10
                grid2(i,j) = snicar_alg10_BBA;
            elseif grid(i,j) >10
                grid2(i,j) = snicar_alg10_BBA;                   
            end
        end
    end
    
    % Plot
    
    BBAlist(counter) = mean2(grid2);
    surf(grid2);
    view(2);
    xlim([1,100])
    ylim([1,100])
    shading interp;
    cbar = colorbar;
    cbar.Limits=[0,1];
    cbar.Label.String = 'Albedo'
    caxis([0,1]);
    drawnow
    %saveas(figure,sprintf('FIG%d.png'));
    
end
