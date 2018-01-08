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
% using a grid-mean. The albedo is used to calculate the radiative forcing
% following Ganey et al (2017: Nat Geoscience) and Dial et al (2018: FEMS)
% by multipling the incoming irradiance by 1-albedo for clean and 'mixed'
% surfaces, the difference between which provides the RF caused by
% impurities.

% One timestep is dimensionless, except that successive SNICAR instances
% have 2x algal biomass in upper layer, meaning 1 timestep is equal to the
% doubling time of the algae, making the time dimension tractable. For a
% doubling time of 3 days, one timestep = 3 days. Albedo change / doubling
% time (days) = albedo change per day. In the final plot, x_time is
% calculated as timestep * doubling_time to enable plotting in units of
% days.

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
       
    % 1. Deal with edge effects (-1s and +1s can push the indexing out of grid
    % range, currently counter starts at 2 and ends at i,j -1, but this
    % leads to a border that doesn't update) and simply isn't plotted.
    
    % 2. Could add an additional randi() dice roll to randomly reduce value
    % of cells, simulating rainfall washaway. Alternatively, a sudden dust
    % deposition event...? These don't have to be random, they could start
    % on a particular ticker value (e.g. at 20th timestep: rain).
        
    % 3. Can the initial grid be set using UAV por satellite imagery?
    
    % 4. How should I deal with dust? constant background? 

    
%%%  VERSION 4 EDITS:
% 1) Updated colormap to show clean snow as white, then grades of red

% 2) Radiative forcing calculations: ***** add details here *****

% 3) Some minor edits to variables names etc for readability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% open a figure for plotting the final albedo grid
figure(1);
figure(2);

% set up empty lists

BBAlist = []; % list to append mean grid value to for albedo grid
spectral_average = [];
IRF_spectral = [];
IRF_sum = [];

%%% USER DEFINED VARIABLES: SET HERE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read incoming irradiance file and assign to variable 'incoming'
fileID = fopen('media/joe/9446-2970/FromGardner/snicar_package/mlw_sfc_flx_frc_clr.txt');
incoming=textscan(fileID,'%f','Delimiter',',');
incoming=cell2mat(incoming); %convert to column vector
fclose(fileID);

% set up time ticker (time_tot = total no of timesteps in run)
time_tot = 25;
timestep = 1;

% %% Set up grid
gridsize = 10000; % total grid area (no. cells)
gridx = 100; % length of x-axis
gridy= 100; % length of y-axis
alg_frac = 0.5; % percentage of algal coverage at start of experiment (all initialise as light algae: class 1)
non_alg_frac = gridsize-alg_frac; % residual = non-algal, assumed clean
doubling_time = 3; % algal doubling time in days (3 is fast, 7 is slow - from literature e.g. Yallop, Stibal) 
chance_insitu = 60; % probability (%) that growth occurs in situ (100 - chance_insitu = chance spreading)

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
      
            if grid(i,j) == 0   % if cell = 0, stay clean: i.e. bloom can't spontaneously form in areas of clean snow
                grid(i,j) = grid(i,j);
                       
            elseif grid(i,j) > 0 && grid(i,j) < 10 % if cell is between 1 and 9, chance of increasing in situ = 'chance_insitu'
                if rand < chance_insitu
                    grid(i,j) = grid(i,j)+1;
                    
                elseif rand > 100-chance_insitu || grid(i,j) == 10 % if rand > 100-chance_insitu OR cell value = 10, bloom spreads, not darkening in situ
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% IRF CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = sum(grid(:)==0) * snicar_clean_spectral; %sum of zeros in grid (weight) * spectral albedo
    B = sum(grid(:)==1) * snicar_alg1_spectral;
    C = sum(grid(:)==2) * snicar_alg2_spectral;
    D = sum(grid(:)==3) * snicar_alg3_spectral;
    E = sum(grid(:)==4) * snicar_alg4_spectral;
    F = sum(grid(:)==5) * snicar_alg5_spectral;
    G = sum(grid(:)==6) * snicar_alg6_spectral;
    H = sum(grid(:)==7) * snicar_alg7_spectral;
    I = sum(grid(:)==8) * snicar_alg8_spectral;
    J = sum(grid(:)==9) * snicar_alg9_spectral;
    K = sum(grid(:)==10) * snicar_alg10_spectral;
    L = sum(grid(:)>10) * snicar_alg10_spectral;
    
    sumspec = A+B+C+D+E+F+G+H+I+J+K+L;
    
% Weighted average spectral albedo    
    spectral_average = sumspec/gridsize;      %% Think really hard about what this really means, i.e. timescale, units
    
    % IRF calculation by multiplying incoming spectral irradiance with area
    % averaged spectral albedo. **Instantaneous** radiative forcing, i,e.
    % NOT integrated over time.
    alb_inv = 1.-spectral_average;
    alb_inv_clean = 1.-snicar_clean_spectral;
    IRF_spectral = alb_inv.*incoming;  %% RF per day for whole area with mixed coverage, incoming = Wm-2, albedo = proportion, space = dimensionless, time = dimensionless (unless assigned)
    IRF_clean = alb_inv_clean.*incoming; % RF per day for whole area but just clean snow
    
    IRF_BBA(counter) = sum(IRF_spectral); % broadband IRF per day (i.e. if algal growth rate = 3d, each timestep = 3d)
    IRF_BBA_clean(counter) = sum(IRF_clean); % braodband IRF per day for clean snow (entire area)

% Plot
    % Set colormap to white-red
    map = [155,0,0
    255, 0, 0
    255, 20, 20
    255, 40, 40
    255, 60, 60
    255, 80, 80
    255, 100, 100
    255, 120, 120
    255, 140, 140
    255, 160, 160
    255, 180, 180
    255, 200, 200
    255, 220, 220
    255, 240, 240
    255, 255, 255];
    map = map./255;
   
    figure(1)
    BBAlist(counter) = mean2(grid2);
    surf(grid2);
    view(2);
    colormap(map);
    xlim([1,100]);
    ylim([1,100]);
    shading interp;
    cbar = colorbar;
    cbar.Limits=[0,0.7];
    cbar.Label.String = 'Albedo'
    caxis([0,0.7]);
    drawnow
    %saveas(figure,sprintf('FIG%d.png'));
    

end
    
    x_time = [1:1:time_tot].*doubling_time; % set up real time in days for plotting
    
    figure(2)
    hold on
    plot(x_time,IRF_BBA,'color','r')
    plot(x_time,BBAlist,'color','b')
    xlabel('Time (days)')
    legend('IRF (W/grid)','Broadband Albedo')
    
    %%% IRF calculations integrated over wavelength and time
   
    RF_BBA = sum(IRF_BBA.*24*60*60); % total RF over area over entire time of model run
    RF_BBA_clean = sum(IRF_BBA_clean.*24*60*60); % total RF over area over entire time of model run
    RF_BBA_impurities = RF_BBA - RF_BBA_clean; % total RF due to impurities over entire area over entire area and time of model run
    
    % daily RF (average daily RF - not calculated for individual days, but calculated by dividing total by no. of days)
    RF_daily = RF_BBA/doubling_time; % RF for whole area per day
    RF_daily_clean = RF_BBA_clean/ doubling_time; % RF for whole area per day, clean snow
    RF_daily_impurities = RF_BBA_impurities / doubling_time; % broadband RF for whole area due to impurities, per day
    
    % daily, per metre (again, calculated as average over space and time
    % divided by gridsize and number of days)
    
    RF_daily_per_m = RF_daily / gridsize; % assuming 1 cell is 1 x 1 m 
    RF_daily_clean_per_m = RF_daily_clean / gridsize; % assuming 1 cell is 1 x 1m
    RF_daily_impurities_per_m = RF_daily_impurities / gridsize; % assuming 1 cell is 1 x 1 m