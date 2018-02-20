% CASPA DRIVER: Joseph Cook, University of Sheffield, UK, Feb 2018

% This code is the master driver script for CASPA including definition of
% initial conditions, setup and model running. This is now the only script
% the user needs to intercat with in order to set up and run CASPA except
% significant modifications are to be made.


%%%%%%%%%%%%%%%%% SECTION 1: CASPA SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here user defined variable values are set and the setup code run.

% Snicar params

GRAIN_SIZE = 0;
BND_TYP = 1;
DIRECT = 1;
APRX_TYP = 1;
DELTA = 1;
coszen = 0.5;
R_sfc = 0.15;

% Initial conditions

dz_init = [0.05,0.05,0.05,0.05,0.05];
rho_snw_init = [200, 200, 200, 200, 300];
rds_snw_init = [350, 350, 350, 350, 350];
rds_coated_init = ["0","0","0","0","0"];
x_init = 0.01e6;
initial_T = [273, 273, 272.75, 272, 270];
fliqs_init = [0,0,0,0,0];
f_refs_init = [0,0,0,0,0];
rho_snw = rho_snw_init;


% run setup code with conditions defined above

CASPA_setup(GRAIN_SIZE, BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc,...
dz_init, rho_snw_init, rds_snw_init, rds_coated_init,...
x_init, initial_T, fliqs_init, f_refs_init);


%%%%%%%%%%%%%%%%%% SECTION 2: CASPA RUN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run CASPA

no_runs = 1;

for i = 1:1:no_runs
    
    [spectral_average, BBAlist, x_time] = CASPA_V01();
    
    BBA(:,i) = BBAlist;

    figure(2)
    plot(x_time,BBAlist,'--k')
    ylabel('Albedo')
    xlabel('Time (days)')
    hold on

end

% plot
if no_runs > 1
    BBA_av = mean(BBA,2)
    figure(2)
    plot(x_time,BBA_av,'color','r','LineWidth',2)
end