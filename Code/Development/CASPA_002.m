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

    % URGENT: Create new SNICAR instances for algal_2 - algal_10
    
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
    
    % 4. See if there is anything you can do to speed this code up - at the
    % moment a grid with 10,000 cells is multiple hours on good
    % laptop (32G RAM, intel i-7 2.8GHz processor)
    
    % 5. Can the initial grid be set using UAV por satellite imagery?
    
    % 6. How should I deal with dust? constant background? 
    
    % 7. Possibility of calculating radiative forcing in each SNICAR run,
    % then calculating energy abalnce and volume of melt,. and linking that
    % to change in grain size? Iterate through and incorporate indirect
    % albedo evolution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a figure for plotting the final albedo grid
figure
hold on

% %% Set up grid
BBAlist = []; % list to append mean grid value to for albedo grid

gridsize = 900; % total grid area (no. cells)
gridx = 30; % length of x-axis
gridy= 30; % length of y-axis
alg_frac = 1; % percentage of algal coverage at start of experiment (all initialise as light algae: class 1)
non_alg_frac = gridsize-alg_frac; % residual = non-algal, assumed clean

% create grid by randomly spacing 1's and 2's in a 10x10 grid
alg_pixels = gridsize*(alg_frac/100) 
non_alg_pixels = gridsize-alg_pixels

grid = reshape([zeros(non_alg_pixels,1) ; ones(alg_pixels,1)],gridx,gridy) ;
grid(:) = grid(randperm(numel(grid))) ;

for counter = 1:1:28    
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
%%%%%%% RUN SNICAR PER CELL AND UPDATE NEW GRID WITH ALBEDO VALUES %%%%%%%%
    
% Note that the remainder of the code is still running within the main
% 'counter' loop, so the surface classification, bloom evolution, call to
% SNICAR and plotting albedo all occur per timestep of the main ticker.

    grid2 = grid; % copy surface class grid for populating with albedos
    
    for i = 1:1:gridx  % Loop through all cells
        for j = 1:1:gridy
            if grid(i,j) == 0  % if class = 0, run clean ice snicar, return BBA and spectral albedo
                [alb_slr,albedo_clean] = SNICAR_driver_clean_CA; 
                grid2(i,j) = alb_slr; % add BBA to appropriate cell in grid2
            elseif grid2(i,j) ==1 % if cell is not 0, run appropriate 'algal' version of snicar
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_1;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==2
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_2;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==3
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_3;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==4
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_4;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==5
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_5;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==6
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_6;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==7
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_7;
                grid2(i,j) = alb_slr;                
            elseif grid(i,j) ==8
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_8;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==9
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_9;
                grid2(i,j) = alb_slr;
            elseif grid(i,j) ==10
                [alb_slr,albedo_alg] = SNICAR_driver_algae_CA_10;
                grid2(i,j) = alb_slr;    
            
            end
        end
    end
    
    BBAlist(counter) = mean2(grid2);
    surf(grid2);
    view(2);
    shading interp;
    cbar = colorbar;
    cbar.Limits=[0,1];
    cbar.Label.String = 'Albedo'
    caxis([0,1]);
    drawnow
    saveas(figure,sprintf('FIG%d.png'));
    pause(0.2)
    
end
