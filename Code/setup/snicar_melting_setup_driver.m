% Written by J Cook (Univ Sheffield) Jan 2018
% This script is is part of the CASPA setup code. Here, initial conditions 
% for the snowpack are set. A modified form of SNICAR with a built-in grain
% evolution parameterisation (snicar_melting_setup.m and
% grain_size_evolution.m) are run, returning updated parameter values.
% These updated values are then used to drive the next SNICAR run, which
% returns updated input values for the next SNICAR run, and so on. 

% Each run, the input values are saved and at the end of the script the
% inpout values used in each successive run is output as one csv file.
% These values can then easil be used to create SNICAR instances in the
% working directory for CASPA (surface classes have a speciifc set of
% SNICAR variable values that are called depending upon the values in each
% cell).
clear

% Set Initial Conditions (at t=0)

% constants

BND_TYP = 1;
DIRECT = 1;
APRX_TYP = 1;
DELTA = 1;
coszen = 0.5;
R_sfc = 0.15;

% variables

dz_init = [0.05,0.05,0.05,0.05,0.05];
rho_snw_init = [100, 200, 200, 300, 300];
rds_snw_init = [100, 150, 200, 250, 300];
rds_coated_init = ["0","0","0","0","0"];
x_init = 0.01e6;
initial_T = [260, 260, 255, 250, 260];
fliqs_init = [0,0,0,0,0];
f_refs_init = [0,0,0,0,0];


% Set up empty variables to append results to

% x_list = zeros(1,1);
% f_refs_list = zeros(1,1);
% T_list = zeros(1,length(f_refs_init));
% rds_list = zeros(1,length(rds_snw_init));


x_list(1,1) = 0;
T_list(1,[1:length(initial_T)]) = initial_T;
f_refs_list(1,[1:length(f_refs_init)]) = f_refs_init;
fliqs_list(1,[1:length(fliqs_init)]) = fliqs_init;
rds_list(1,[1:length(rds_snw_init)]) = rds_snw_init;



%%%%%%% BEGIN RUN 1 %%%%%%%%%%%%%
[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);

% add new variable values to lists. Each run is a separate row, each snow layer
% is a separate column
x_list(2,1) = x_init;
T_list(2,[1:length(new_T)]) = new_T;
f_refs_list(2,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(2,[1:length(fliqs_init)]) = new_fliqs;
rds_list(2,[1:length(rds_snw_init)]) = round(new_r);

% update outputs for next run
x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%% BEGIN RUN 2  %%%%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);


x_list(3,1) = x_init;
T_list(3,[1:length(new_T)]) = new_T;
f_refs_list(3,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(3,[1:length(fliqs_init)]) = new_fliqs;
rds_list(3,[1:length(rds_snw_init)]) = round(new_r);

Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);
x_init = x_init*2;


%%%%%%%% BEGIN RUN 3 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(4,1) = x_init;
T_list(4,[1:length(new_T)]) = new_T;
f_refs_list(4,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(4,[1:length(fliqs_init)]) = new_fliqs;
rds_list(4,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 4 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(5,1) = x_init;
T_list(5,[1:length(new_T)]) = new_T;
f_refs_list(5,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(5,[1:length(fliqs_init)]) = new_fliqs;
rds_list(5,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 5 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(6,1) = x_init;
T_list(6,[1:length(new_T)]) = new_T;
f_refs_list(6,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(6,[1:length(fliqs_init)]) = new_fliqs;
rds_list(6,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 6 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(7,1) = x_init;
T_list(7,[1:length(new_T)]) = new_T;
f_refs_list(7,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(7,[1:length(fliqs_init)]) = new_fliqs;
rds_list(7,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 7 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(8,1) = x_init;
T_list(8,[1:length(new_T)]) = new_T;
f_refs_list(8,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(8,[1:length(fliqs_init)]) = new_fliqs;
rds_list(8,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 8 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(9,1) = x_init;
T_list(9,[1:length(new_T)]) = new_T;
f_refs_list(9,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(9,[1:length(fliqs_init)]) = new_fliqs;
rds_list(9,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 9 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(10,1) = x_init;
T_list(10,[1:length(new_T)]) = new_T;
f_refs_list(10,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(10,[1:length(fliqs_init)]) = new_fliqs;
rds_list(10,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%% BEGIN RUN 10 %%%%%%%%%%

[new_r, new_f_refs, new_T, new_fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
x_list(11,1) = x_init;
T_list(11,[1:length(new_T)]) = new_T;
f_refs_list(11,[1:length(f_refs_init)]) = new_f_refs;
fliqs_list(11,[1:length(fliqs_init)]) = new_fliqs;
rds_list(11,[1:length(rds_snw_init)]) = round(new_r);

x_init = x_init*2;
Initial_T = new_T;
f_refs_init = new_f_refs;
fliqs_init = new_fliqs;
rds_snw_init = round(new_r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CREATE OUTPUT CSV %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create labelled table for each output variable

T_table = array2table(T_list,'VariableNames',{'T_top','T_z2','T_z3','T_z4','T_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
x_table = array2table(x_list,'VariableNames',{'x'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'})
f_ref_table = array2table(f_refs_list,'VariableNames',{'f_ref_top','f_ref_z2','f_ref_z3','f_ref_z4','f_ref_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
fliqs_table = array2table(fliqs_list,'VariableNames',{'fliqs_top','fliqs_z2','fliqs_z3','fliqs_z4','fliqs_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
rds_table = array2table(rds_list,'VariableNames',{'r_top','r_z2','r_z3','r_z4','r_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});

%concatenate tables into one large output table

output_table = [x_table,rds_table,fliqs_table,f_ref_table,T_table]


% Write table to csv

%writetable(output_table,'CASPA_SNICAR_Inputs.csv','WriteRowNames',true) 