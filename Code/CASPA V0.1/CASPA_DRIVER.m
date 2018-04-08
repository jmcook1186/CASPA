% CASPA DRIVER: Joseph Cook, University of Sheffield, UK, Feb 2018

% This code is the master driver script for CASPA including definition of
% initial conditions, setup and model running. This is now the only script
% the user needs to interact with in order to set up and run CASPA except
% significant modifications are to be made.


clear

%%%%%%%%%%%%%%%%% SECTION 1: CASPA SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here user defined variable values are set and the setup code run.

folder_path = "/media/joe/9446-2970/FromGardner/snicar_package/";

% Snicar params

GRAIN_SIZE = 1;   % switch grain evolution model on/off (1/0)
BND_TYP = 1;      % 1 = spectral
DIRECT = 1;       % Direct illumination = 1, diffuse =  0
APRX_TYP = 1;     % Eddington (1), quadrature (2) or hemispheric mean(3)
DELTA = 1;        % 1 = apply delta approximation for anisotropic scattering
coszen = 0.5;     % cosine of solar zenith angle
R_sfc = 0.15;     % albedo of underlying surface


% %% Set up grid

time_tot = 30;  % 30 timesteps = 90 days if doubling time = 3
timestep = 1;
gridsize = 40000; % total grid area (no. cells)
length_scale = 0.1; % define length of each pixel in metres
alg_frac = 0.5; % percentage of algal coverage at start of experiment (all initialise as light algae: class 1)
non_alg_frac = gridsize-alg_frac; % residual = non-algal, assumed clean
doubling_time = 3; % algal doubling time in days (3 is fast, 7 is slow - from literature e.g. Yallop, Stibal) 
chance_insitu = 55; % probability (%) that growth occurs in situ (100 - chance_insitu = chance spreading)

% Initial conditions

x_init = 5000; % algal concentration at t=0
y_init = 500000; % dust concentration at t=0
scavenging_rate = 1.1; % how much scavenging of dust occurs per timestep (2 = double dust per timestep)
dz_init = [0.05,0.05,0.05,0.05,0.05];
rho_snw_init = [150, 150, 150, 150, 150];
rds_snw_init = [100, 200, 200, 200, 200];
rds_coated_init = ["0","0","0","0","0"];
initial_T = [273.15, 273, 272, 270, 268];
fliqs_init = [0,0,0,0,0];
f_refs_init = [0,0,0,0,0];
rho_snw = rho_snw_init;

% run setup code with conditions defined above

CASPA_setup(GRAIN_SIZE, BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init,x_init, y_init, scavenging_rate, initial_T, fliqs_init, f_refs_init,folder_path);


%%%%%%%%%%%%%%%%%% SECTION 2: CASPA RUN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_runs = 5;

for i = 1:1:num_runs
    
    [spectral_average, BBAlist, x_time, Tot_biomass_per_m, Tot_alg_pixels_percent_list] = CASPA_V01(time_tot,timestep, gridsize, length_scale, alg_frac, doubling_time,chance_insitu, folder_path,rho_snw);
    
    BBA(:,i) = BBAlist;
    Biomass(:,i) = Tot_biomass_per_m;

    figure(2)
    yyaxis right
   
    plot(x_time,BBAlist,'--m','LineWidth',0.1)
    ylabel('Albedo')
    xlabel('Time (days)')
    hold on
    
    yyaxis left
    plot(x_time,Tot_biomass_per_m,'--c','LineWidth',0.1)
    ylabel('Biomass g/m^2')

    figure(3)
    yyaxis right
    
    plot(x_time,Tot_alg_pixels_percent_list,'--m','color','r','LineWidth',0.2)
    xlabel('Time (days)')
    ylabel('% coverage')
    
    yyaxis left
    plot(x_time,BBAlist, '--c','LineWidth',0.2)
    ylabel('albedo')
    hold on
    
end

% plot

if num_runs > 1
    
    BBA_av = mean(BBA,2);
    Biomass_av = mean(Biomass,2);
    
    figure(2)
    yyaxis right
    plot(x_time,BBA_av,'color','r','LineWidth',2)
    
    yyaxis left
    plot(x_time, Biomass_av,'color','b','LineWidth',2)
    
   
    
end