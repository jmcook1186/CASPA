% CASPA: CELLULAR AUTOMATON FOR SNOWPACK ALBEDO
% Joseph Cook, University of SHeffield, February 2018

% %% OVERVIEW

% This model simulates a snow surface with initial physical/optical properties
% defined by the user, then predicted using an algal growth and snow grain evolution
% model. This surface is represented as a cellular automaton,
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
% scheme BioSNICAR. These instances are pre-run by the setup code at the
% same time as running the snow grain evolution model, with input
% values automatically updated as the snow and algal bloom develop over
% time. The broadband albedo and spectral albedo redicted by SNICAR are
% saved to .mat files in the working directory. 

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
        
    % 2. Can the initial grid be set using UAV or satellite imagery?
    
    % 3. How should I deal with dust? constant background? 
    
    % 4. Need to reduce IRF etc due to day/night irradiance changes?
    % Currently assumes constant irradiance throughout timestep (no night
    % time)

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
time_tot = 30;
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

% LOAD broadband and spectral albedos generated using SNICAR in the setup
% driver code. This is much faster than calling SNICAR in CASPA and means no
% manual entry of input values is required).

SPEC = load('spectral_list.mat');
BBA = load('BBA.mat');
      
%%% Loop through grid and replace surface class with albedo from SNICAR lookup
% library populated in loop above.
        
    for i = 1:1:gridx  % Loop through all cells
        for j = 1:1:gridy
            if grid(i,j) == 0  % if class = 0, run clean ice snicar, return BBA and spectral albedo 
                grid2(i,j) = BBA.BBA(1); % add BBA to appropriate cell in grid2
            elseif grid2(i,j) ==1 % if cell is not 0, run appropriate 'algal' version of snicar
                grid2(i,j) = BBA.BBA(2);
            elseif grid(i,j) ==2
                grid2(i,j) = BBA.BBA(3);
            elseif grid(i,j) ==3
                grid2(i,j) = BBA.BBA(4);
            elseif grid(i,j) ==4
                grid2(i,j) = BBA.BBA(5);
            elseif grid(i,j) ==5
                grid2(i,j) = BBA.BBA(6);
            elseif grid(i,j) ==6
                grid2(i,j) = BBA.BBA(7);
            elseif grid(i,j) ==7
                grid2(i,j) = BBA.BBA(8);                
            elseif grid(i,j) ==8
                grid2(i,j) = BBA.BBA(9);
            elseif grid(i,j) ==9
                grid2(i,j) = BBA.BBA(10);
            elseif grid(i,j) ==10
                grid2(i,j) = BBA.BBA(11);
            elseif grid(i,j) >10
                grid2(i,j) = BBA.BBA(11);                   
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% IRF CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

    SPEC.spectral_list(:,2)
    
    A = sum(grid(:)==0) * SPEC.spectral_list(:,1); %sum of zeros in grid (weight) * spectral albedo
    B = sum(grid(:)==1) * SPEC.spectral_list(:,2);
    C = sum(grid(:)==2) * SPEC.spectral_list(:,3);
    D = sum(grid(:)==3) * SPEC.spectral_list(:,4);
    E = sum(grid(:)==4) * SPEC.spectral_list(:,5);
    F = sum(grid(:)==5) * SPEC.spectral_list(:,6);
    G = sum(grid(:)==6) * SPEC.spectral_list(:,7);
    H = sum(grid(:)==7) * SPEC.spectral_list(:,8);
    I = sum(grid(:)==8) * SPEC.spectral_list(:,9);
    J = sum(grid(:)==9) * SPEC.spectral_list(:,10);
    K = sum(grid(:)==10) * SPEC.spectral_list(:,11);
    L = sum(grid(:)>10) * SPEC.spectral_list(:,11);
    
    sumspec = A+B+C+D+E+F+G+H+I+J+K+L;
    
% Weighted average spectral albedo    
    spectral_average = sumspec/gridsize;      %% Think really hard about what this really means, i.e. timescale, units
    
    % IRF calculation by multiplying incoming spectral irradiance with area
    % averaged spectral albedo. **Instantaneous** radiative forcing, i,e.
    % NOT integrated over time.
    alb_inv = 1.-spectral_average;
    alb_inv_clean = 1.-SPEC.spectral_list(:,1);
    IRF_spectral = alb_inv.*incoming;  %% RF per day for whole area with mixed coverage, incoming = Wm-2, albedo = proportion, space = dimensionless, time = dimensionless (unless assigned)
    IRF_clean = alb_inv_clean.*incoming; % RF per day for whole area but just clean snow
    
    IRF_BBA(counter) = sum(IRF_spectral); % broadband IRF at each timestep
    IRF_BBA_clean(counter) = sum(IRF_clean); % broadband IRF at each timestep for clean snow (entire area)

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
    cbar.Label.String = 'Albedo';
    caxis([0,0.7]);
    drawnow
    %saveas(figure,sprintf('FIG%d.png'));
    

end
    
    %%% IRF calculations integrated over wavelength and time
    
    % NOTE that the IRF here is measured per second by default. This is
    % because the radiative transfer equations are calculations of flux,
    % where flux is energy per unit time measured in Watts, where 1 watt =
    % 1 joule/second. 
    
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
    
    
    % Plot Fig 2 (RF and Albedo against time)
    
    x_time = [1:1:time_tot].*doubling_time; % set up real time in days for plotting
    
    figure(2)
    hold on
    yyaxis right
    plot(x_time,IRF_BBA,'color','r')
    ylabel('IRF (W/grid)')
    yyaxis left
    plot(x_time,BBAlist,'color','b')
    ylabel('Albedo')
    xlabel('Time (days)')
    legend('IRF (W/grid)','Broadband Albedo')
    