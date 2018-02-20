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

% Note that the new effective radii added to the table are the overall new
% radii which take into account the grain evolution and also the refrozen
% fraction in each layer (using a weighted mean, assuming refrozen water
% has an effective radius of 1500 microns).

function[] = CASPA_setup(GRAIN_SIZE, BND_TYP,DIRECT,APRX_TYP,DELTA,...
coszen,R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init,...
x_init, initial_T, fliqs_init, f_refs_init)

% append initial conditions to 1st line of output table
x_list(1,1) = 0;
T_list(1,[1:length(initial_T)]) = initial_T;
f_refs_list(1,[1:length(f_refs_init)]) = f_refs_init;
fliqs_list(1,[1:length(fliqs_init)]) = fliqs_init;
rds_list(1,[1:length(rds_snw_init)]) = rds_snw_init;
dz_list(1,[1:length(dz_init)]) = dz_init;


if GRAIN_SIZE==1
    %%%%%%% BEGIN RUN 1 %%%%%%%%%%%%%
    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);

    % add new variable values to lists. Each run is a separate row, each snow layer
    % is a separate column
    x_list(2,1) = x_init;
    T_list(2,[1:length(new_T)]) = new_T;
    f_refs_list(2,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(2,[1:length(fliqs_init)]) = fliqs;
    rds_list(2,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(2,[1:length(dz_init)]) = dz_init;

    % update outputs for next run
    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%% BEGIN RUN 2  %%%%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);


    x_list(3,1) = x_init;
    T_list(3,[1:length(new_T)]) = new_T;
    f_refs_list(3,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(3,[1:length(fliqs_init)]) = fliqs;
    rds_list(3,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(3,[1:length(dz_init)]) = dz_init;

    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    x_init = x_init*2;
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 3 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(4,1) = x_init;
    T_list(4,[1:length(new_T)]) = new_T;
    f_refs_list(4,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(4,[1:length(fliqs_init)]) = fliqs;
    rds_list(4,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(4,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 4 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(5,1) = x_init;
    T_list(5,[1:length(new_T)]) = new_T;
    f_refs_list(5,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(5,[1:length(fliqs_init)]) = fliqs;
    rds_list(5,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(5,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 5 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(6,1) = x_init;
    T_list(6,[1:length(new_T)]) = new_T;
    f_refs_list(6,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(6,[1:length(fliqs_init)]) = fliqs;
    rds_list(6,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(6,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 6 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(7,1) = x_init;
    T_list(7,[1:length(new_T)]) = new_T;
    f_refs_list(7,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(7,[1:length(fliqs_init)]) = fliqs;
    rds_list(7,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(7,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 7 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(8,1) = x_init;
    T_list(8,[1:length(new_T)]) = new_T;
    f_refs_list(8,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(8,[1:length(fliqs_init)]) = fliqs;
    rds_list(8,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(8,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 8 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(9,1) = x_init;
    T_list(9,[1:length(new_T)]) = new_T;
    f_refs_list(9,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(9,[1:length(fliqs_init)]) = fliqs;
    rds_list(9,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(9,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 9 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(10,1) = x_init;
    T_list(10,[1:length(new_T)]) = new_T;
    f_refs_list(10,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(10,[1:length(fliqs_init)]) = fliqs;
    rds_list(10,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(10,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%% BEGIN RUN 10 %%%%%%%%%%

    [overall_new_r, f_refs, new_T, fliqs, dz_new] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz_init, rho_snw_init, rds_snw_init, rds_coated_init, x_init, initial_T, fliqs_init, f_refs_init);
    x_list(11,1) = x_init;
    T_list(11,[1:length(new_T)]) = new_T;
    f_refs_list(11,[1:length(f_refs_init)]) = f_refs;
    fliqs_list(11,[1:length(fliqs_init)]) = fliqs;
    rds_list(11,[1:length(rds_snw_init)]) = round(overall_new_r);
    dz_list(11,[1:length(dz_init)]) = dz_init;

    x_init = x_init*2;
    Initial_T = new_T;
    f_refs_init = f_refs;
    fliqs_init = fliqs;
    rds_snw_init = round(overall_new_r);
    dz_init = dz_new;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% CREATE OUTPUT CSV %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create labelled table for each output variable

    T_table = array2table(T_list,'VariableNames',{'T_top','T_z2','T_z3','T_z4','T_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    x_table = array2table(x_list,'VariableNames',{'x'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    f_ref_table = array2table(f_refs_list,'VariableNames',{'f_ref_top','f_ref_z2','f_ref_z3','f_ref_z4','f_ref_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    fliqs_table = array2table(fliqs_list,'VariableNames',{'fliqs_top','fliqs_z2','fliqs_z3','fliqs_z4','fliqs_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    rds_table = array2table(rds_list,'VariableNames',{'r_top','r_z2','r_z3','r_z4','r_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    dz_table = array2table(dz_list,'VariableNames',{'dz_top','dz_z2','dz_z3','dz_z4','dz_btm'},'RowNames',{'SNICAR_INIT','SNICAR_1','SNICAR_2','SNICAR_3','SNICAR_4','SNICAR_5','SNICAR_6','SNICAR_7','SNICAR_8','SNICAR_9','SNICAR_10'});
    %concatenate tables into one large output table (currently only relevant
    %variables for creating snicar instances (x, rds, water coat thickness)

    output_table = [x_table,T_table,dz_table,rds_table]

    snicar_params = [BND_TYP,DIRECT,APRX_TYP,DELTA,coszen,R_sfc];
    snicar_params_table = array2table(snicar_params,'VariableNames',{'BND_TYP','DIRECT','APRX_TYP','DELTA','coszen','R_sfc'})


    % Write tables to csv

    %writetable(output_table,'CASPA_SNICAR_Inputs.csv','WriteRowNames',true) 
    %writetable[snicar_params_table,'CASPA,SNICAR_params.csv']

    % Uncomment the call to mie_driver_caspa in the following loop to run the
    % mie solver and add the netcdf files for the relevant new grain sizes 
    % and water coatings to the working directory. Pads with zeros if the no of
    % digits is <4, meaning the filenaming convention is consistent with the
    % SNICAR driver in matlab and FORTRAN and also CLM.


    for i = 1:1:numel(rds_list)

        rad = rds_list(i);

        % set path

        folder_path = '/media/joe/9446-2970/FromGardner/snicar_package/';

        a = sprintf('%04s',num2str(rad));

        % if water coating required, add 'coated' into filename 

        filename = strcat('ice_wrn_',a,'.nc');
        fullpath = strcat(folder_path,filename);

        % check to see if file is already in the working directory, if it is
        % just raise a flag and do not run the mie solver. If the file does not
        % already exist, run the mie solver and add file to wdir

        if exist(fullpath,'file')

            A = strcat('NetCDF for ', filename, ' already exists');

        else

            A = sprintf('doing mie calc for %s',filename)
            mie_driver_caspa(rad);

        end

    end


    % Run SNICAR with input values from the setup code above. The
    % 'output_table' from the setup code provides the inout values for each
    % instance of SNICAR. Two files are saved to the workspace: spectral_list
    % and BBA. Each column of 'spectral_list.mat' contains the spectral albedo
    % of the snow for each SNICAR instance (1-11). The file 'BBA.mat' contains
    % the broadband albedo calculated by SNICAR for each instance. These can
    % then be read by CASPA instead of running SNICAR each time. No manual
    % inputting of values into SNICAR is required using this method = major
    % time saving and avoids human errors.

    spectral_list = zeros(470,length(x_list));

    for i = 1:1:length(x_list)
        dz = dz_list(i,:);
        rho_snw(i,:) = rho_snw_init;
        rds_snw(i,:) = rds_list(i,:);
        x = x_list(i,:);

        [BBalb,spectral] = SNICAR_function_CASPA(BND_TYP,DIRECT,APRX_TYP,DELTA,coszen,R_sfc,dz,rho_snw,rds_snw,x);

        BBA(i) = BBalb;
        spectral_list(:,i) = spectral;

    end
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% IF GRAIN SIZE EVOLUTION MODEL SWITCHED OFF %%%%%%%%%%%%%%%%%%%%


else
    spectral_list = zeros(470,11);
    x_list(1,1) = 0;
    x_list(1,2) = x_init;
    
    for i = 1:1:11

        dz_list(i,[1:length(dz_init)]) = dz_init;
        rho_snw(i,:) = rho_snw_init;
        rds_list(i,[1:length(rds_snw_init)]) = rds_snw_init;
        
        if i > 1
            x_list(1,i+1) = x_list(i)*2;
        end
        
    end
    
    
    for i = 1:1:11
        
        x = x_list(1,i);
        dz = dz_init;
        rho_snw = rho_snw_init;
        rds_snw = rds_snw_init;
        
        [BBalb,spectral] = SNICAR_function_CASPA(BND_TYP,DIRECT,APRX_TYP,DELTA,coszen,R_sfc,dz,rho_snw,rds_snw,x);

        BBA(i) = BBalb;
        spectral_list(:,i) = spectral;
    
    end

x_list = x_list(1:11);

end

    
save('spectral_list.mat','spectral_list');
save('BBA.mat','BBA');
save('X_list.mat','x_list')
save('dz.mat','dz_list')




end
