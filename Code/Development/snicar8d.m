% The Snow, Ice, and Aerosol Radiative (SNICAR) Model: written by Mark
% Flanner (flanner@umich.edu)
%
% This function calculates radiative fluxes at the interfaces of
% multiple snow layers using a 2-stream approximation with a
% tri-diagonal matrix solution, presented by Toon, O.B. et al.,
% Journal for Geophysical Research, vol 94 D13 p 16287, 1989.  Thus, a
% vertically inhomogenous medium can be represented.
%
% The function reads in Mie optical properties from external NetCDF
% files based on defined properties of each layer.
%
% Example input parameters are defined below. The user has option of
% running the program as a script. To do so, comment out the function
% line and set '1==1' below:
%
%
%%%%%%%%%%  Input parameters: %%%%%%%%%%%
% BND_TYP:      Spectral grid (=1 for 470 bands. This is the
%               only functional option in this distribution)
% DIRECT:       Direct or diffuse incident radiation (1=direct, 0=diffuse)
% APRX_TYP:     Two-Stream Approximation Type:
%                  1 = Eddington
%                  2 = Quadrature
%                  3 = Hemispheric Mean
%                  NOTE: Delta-Eddington Approximation is probably
%                  best for the visible spectrum, but can provide
%                  negative albedo in the near-IR under diffuse
%                  light conditions. Hence, the Hemispheric Mean
%                  approximation is recommended for general use.
% DELTA:        1=Use Delta approximation (Joseph, 1976), 0=don't use
% coszen:       cosine of solar zenith angle (only applies when DIRECT=1)
% R_sfc:        broadband albedo of underlying surface 
%                  (user can also define a spectrally-dependent albedo below)
% dz:           array of snow layer thicknesses [meters]. Top layer has index 1. 
%                  The length of this array defines number of snow layers
% rho_snw:      array of snow layer densities [kg m-3]. Must have same length as dz
% rds_snw:      array of snow layer effective grain radii [microns]. Must have same length as dz
% nbr_aer:      number of aerosol species in snowpack
% mss_cnc_sot1: mass mixing ratio of black carbon species 1 (uncoated BC)
%                 (units of parts per billion, ng g-1)
% mss_cnc_sot2: mass mixing ratio of black carbon species 2 (sulfate-coated BC)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst1: mass mixing ratio of dust species 1 (radii of 0.05-0.5um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst2: mass mixing ratio of dust species 2 (radii of 0.5-1.25um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst3: mass mixing ratio of dust species 3 (radii of 1.25-2.5um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst4: mass mixing ratio of dust species 4 (radii of 2.5-5.0um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_ash1: mass mixing ratio of volcanic ash species 1
%                 (units of parts per billion, ng g-1)
% fl_sot1:      name of file containing optical properties for BC species 1
% fl_sot2:      name of file containing optical properties for BC species 2
% fl_dst1:      name of file containing optical properties for dust species 1
% fl_dst2:      name of file containing optical properties for dust species 2
% fl_dst3:      name of file containing optical properties for dust species 3
% fl_dst4:      name of file containing optical properties for dust species 4
% fl_ash1:      name of file containing optical properties for ash species 1
%
% NOTE:  The spectral distribution of incident radiation is, by
% default, one typical of a mid-latitude winter. The user can
% change the incident flux distribution by changing the text files
% loaded below under "Set incident solar flux spectral distribution:"
%
%

%%%%%  Output data: %%%%%

% All output data is contained in the matrix "data_out"
% data_out(:,1) = wavelength (microns, m^-6)
% data_out(:,2) = spectrally-dependent albedo;
% data_out(1,3) = broadband (0.3-5.0um) solar albedo
% data_out(2,3) = visible (0.3-0.7um) albedo;
% data_out(3,3) = near-IR (0.7-5.0um) albedo;
% data_out(4,3) = total solar absorption by snowpack (W/m2)
%                 (not including underlying substrate)
% data_out(5,3) = solar absorption in top snow layer (W/m2)
% data_out(6,3) = solar absorption in second snow layer (W/m2)
% (see "data_out" fields near end of script for more output and to
% add new output)




function data_out = snicar8d(BND_TYP_IN, DIRECT_IN, APRX_TYP_IN, ...
                             DELTA_IN, coszen_in, R_sfc_in, dz_in, ...
                             rho_snw_in, rds_snw_in, rds_coated_in, nbr_aer_in, ...
                             mss_cnc_sot1_in, mss_cnc_sot2_in, ...
                             mss_cnc_dst1_in, mss_cnc_dst2_in, ...
                             mss_cnc_dst3_in, mss_cnc_dst4_in, ...
                             mss_cnc_ash1_in,mss_cnc_bio1_in, mss_cnc_bio2_in,mss_cnc_bio3_in,mss_cnc_bio4_in,mss_cnc_bio5_in, mss_cnc_bio6_in, mss_cnc_bio7_in, mss_cnc_RBio1_in, mss_cnc_hematite_in, mss_cnc_mixed_sand_in, fl_sot1_in, fl_sot2_in, ...
                             fl_dst1_in, fl_dst2_in, fl_dst3_in, ...
                             fl_dst4_in, fl_ash1_in,fl_bio1_in,fl_bio2_in,fl_bio3_in,fl_bio4_in, fl_bio5_in,fl_bio6_in, fl_bio7_in, fl_RBio1_in, fl_hematite_in, fl_mixed_sand_in);

    
if (1==0)
    % DEFINE ALL INPUT HERE (not called from external function);
    
    clear;

    % SNOW LAYER THICKNESSES [m]:
    dz = [0.02 9.98];
 
    nbr_lyr = length(dz);  % number of snow layers
  
  
    % SNOW DENSITY FOR EACH LAYER (units: kg/m3)
    rho_snw(1:nbr_lyr) = 150;  
  

    % SNOW GRAIN SIZE FOR EACH LAYER (units: microns):
    rds_snw(1:nbr_lyr) = 100;
  

    % ALTERNATIVE SNOW GRAIN SIZE: 9999
    % SET OPTICAL PROPERTIES FROM USER-DEFINED FILE:
    %rds_snw(1:nbr_lyr)=9999;
    %fl_in_usr = 'H:FromGardner/snicar_package/ice_wrn_20000.nc';

    % ALTERNATIVE GRAIN SIZES FOR COATED SPHERES
    
    rds_coated(1:nbr_lyr) = 100;   % JC EDIT
    
    % NUMBER OF PARTICLE SPECIES IN SNOW (ICE EXCLUDED)
    
    nbr_aer = 17; % JC EDIT
  
    
    % PARTICLE MASS MIXING RATIOS (units: ng g-1)
    mss_cnc_sot1(1:nbr_lyr)  = 100.0;  % uncoated black carbon
    mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % coated black carbon
    mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust species 1
    mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust species 2
    mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust species 3
    mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust species 4
    mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash species 1
    mss_cnc_bio1(1:nbr_lyr)  = 0.0;    % biological impurity 1 % JC EDIT
    mss_cnc_bio2(1:nbr_lyr)  = 0.0;    % biological impurity 2 % JC EDIT
    mss_cnc_bio3(1:nbr_lyr)  = 0.0;    % biological impurity 3 % JC EDIT
    mss_cnc_bio4(1:nbr_lyr)  = 0.0;    % biological impurity 4 % JC EDIT
    mss_cnc_bio5(1:nbr_lyr)  = 0.0;    % biological impurity 5 % JC EDIT
    mss_cnc_bio6(1:nbr_lyr)  = 0.0;    % biological impurity 6 % JC EDIT
    mss_cnc_bio7(1:nbr_lyr)  = 0.0;    % biological impurity 7 % JC EDIT
    mss_cnc_RBio1(1:nbr_lyr) = 0.0; % Algae measured pigments % JC EDIT
    mss_cnc_hematite(1:nbr_lyr)= 0.0;    % hematite type 1 % JC EDIT
    mss_cnc_mixed_sand(1:nbr_lyr) = 0.0; % mixed sand % JC EDIT
    
    % FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
    % (ideally, these files should exist in all 'band' directories)
    fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
    fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
    fl_dst1  = 'aer_dst_bln_20060904_01.nc';
    fl_dst2  = 'aer_dst_bln_20060904_02.nc';
    fl_dst3  = 'aer_dst_bln_20060904_03.nc';
    fl_dst4  = 'aer_dst_bln_20060904_04.nc';
    fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
    fl_bio1  = 'biological_1.nc' % Biological impurity 1 (30um diameter, 1.5%chll a, 10% each 1 % 2 carotenoids) % JC EDIT
    fl_bio2  = 'biological_2.nc'; % Biological impurity 2 (30um diameter, 1.5%chll a, 5% each 1 % 2 carotenoids) % JC EDIT
    fl_bio3  = 'biological_3.nc'; % Biological impurity 3 (30um diameter, 1.5%chll a, 1% each 1 % 2 carotenoids) % JC EDIT
    fl_bio4  = 'biological_4.nc'; % Biological impurity 4 (30um diameter, 1.5%chll a only) % JC EDIT
    fl_bio5  = 'biological_5.nc'; % Biological impurity 5 (10um diameter, pigs as per bio2) % JC EDIT
    fl_bio6  = 'biological_6.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2) % JC EDIT
    fl_bio7  = 'biological_7.nc'; % Biological impurity 6 (20um diameter, pigs as per bio2) % JC EDIT
    fl_RBio1 = 'Real_Bio1.nc' % Biological impurity with measured pigments, 20 micron diam % JC EDIT
    fl_hematite = 'Hematite.nc'; % hematite type 1 % JC EDIT
    fl_mixed_sand = 'mixed_sand.nc'; % mixed sand % JC EDIT
    
    % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
    
    coszen=0.5;

  
    % REFLECTANCE OF SURFACE UNDERLYING SNOW:
    % (value applied to all wavelengths.  user can also specify
    % spectrally-dependent ground albedo below)
    
    R_sfc_all_wvl = 0.25;
  
  
    % RADIATIVE TRANSFER CONFIGURATION:
    BND_TYP  = 1;        % 1= 470 spectral bands, 2= 5 bands, 3= 3 bands
    DIRECT   = 1;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 3;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta
  
    
   
else
    % routine has been called from external function, with needed
    % input paramaters
    
    % RADIATIVE TRANSFER CONFIGURATION:
    BND_TYP       = BND_TYP_IN;   % 1= 470 spectral bands, 2= 5 bands, 3= 3 bands
    DIRECT        = DIRECT_IN;    % 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP      = APRX_TYP_IN;  % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA         = DELTA_IN;     % 1= Apply Delta approximation, 0= No delta
  
    coszen        = coszen_in;
    R_sfc_all_wvl = R_sfc_in;

    dz            = dz_in;
    nbr_lyr       = length(dz);
    rho_snw       = rho_snw_in;
    rds_snw       = rds_snw_in;
    nbr_aer       = nbr_aer_in;
    rds_coated    = rds_coated_in;

    mss_cnc_sot1  = mss_cnc_sot1_in;
    mss_cnc_sot2  = mss_cnc_sot2_in;  
    mss_cnc_dst1  = mss_cnc_dst1_in;  
    mss_cnc_dst2  = mss_cnc_dst2_in;  
    mss_cnc_dst3  = mss_cnc_dst3_in;  
    mss_cnc_dst4  = mss_cnc_dst4_in;  
    mss_cnc_ash1  = mss_cnc_ash1_in;  
    mss_cnc_bio1  = mss_cnc_bio1_in;   % JC EDIT
    mss_cnc_bio2  = mss_cnc_bio2_in;   % JC EDIT
    mss_cnc_bio3  = mss_cnc_bio3_in;   % JC EDIT
    mss_cnc_bio4  = mss_cnc_bio4_in;   % JC EDIT
    mss_cnc_bio5  = mss_cnc_bio5_in;   % JC EDIT
    mss_cnc_bio6  = mss_cnc_bio6_in;   % JC EDIT
    mss_cnc_bio7  = mss_cnc_bio7_in;   % JC EDIT
    mss_cnc_RBio1 = mss_cnc_RBio1_in; % JC EDIT
    mss_cnc_hematite = mss_cnc_hematite_in;   % JC EDIT
    mss_cnc_mixed_sand = mss_cnc_mixed_sand_in; % JC EDIT
    
    fl_sot1       = fl_sot1_in;
    fl_sot2       = fl_sot2_in;
    fl_dst1       = fl_dst1_in;
    fl_dst2       = fl_dst2_in;
    fl_dst3       = fl_dst3_in;
    fl_dst4       = fl_dst4_in;
    fl_ash1       = fl_ash1_in;
    fl_bio1       = fl_bio1_in;   % JC EDIT
    fl_bio2       = fl_bio2_in;   % JC EDIT
    fl_bio3       = fl_bio3_in;   % JC EDIT
    fl_bio4       = fl_bio4_in;   % JC EDIT
    fl_bio5       = fl_bio5_in;   % JC EDIT
    fl_bio6       = fl_bio6_in;   % JC EDIT
    fl_bio7       = fl_bio7_in;   % JC EDIT
    fl_RBio1      = fl_RBio1_in; % JC EDIT
    fl_hematite   = fl_hematite_in;   % JC EDIT
    fl_mixed_sand = fl_mixed_sand_in; % JC EDIT

    
end;


% Set snow Mie directory based on band number (and direct or
% diffuse flux):
if (BND_TYP==1)
    wrkdir='./';
elseif(BND_TYP==2)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    if (DIRECT==1)
        wrkdir='/data/flanner/mie_clm/mie_clm_drc/';
    elseif (DIRECT==0)
        wrkdir='/data/flanner/mie_clm/mie_clm_dfs/';
    end;
elseif(BND_TYP==3)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    if (DIRECT==1)
        wrkdir='/data/flanner/mie_clm/mie_clm_drc_3bnd/';
    elseif (DIRECT==0)
        wrkdir='/data/flanner/mie_clm/mie_clm_dfs_3bnd/';
    end;
end;

% Set wavelength grid (um) and wavelength number:
wvl     = ncread(strcat(wrkdir,'ice_wrn_0100.nc'),'wvl').*10^6;
nbr_wvl = length(wvl);


% REFLECTANCE OF UNDERLYING SURFACE
% (optional: user can set spectrally-dependent surface albedo here):
R_sfc(1:nbr_wvl,1) = R_sfc_all_wvl;


% file substrings for snow/aerosol Mie parameters:
% Warren and Brandt (2008) ice optical properties:
fl_stb1 = 'ice_wrn_';
fl_stb2 = '.nc';
fl_stb3 = 'coated_';

% Calculate mu_not
%mu_not=cos((slr_znt/360)*2*pi);
mu_not = coszen;


% Set incident solar flux spectral distribution:
if (BND_TYP==1)

    if (DIRECT == 1)
        % mid-latitude winter, clear-sky:
        load mlw_sfc_flx_frc_clr.txt;
        mlw_sfc_flx_frc_clr(find(mlw_sfc_flx_frc_clr==0))=1e-30;
        flx_slr = mlw_sfc_flx_frc_clr;
        
        % Summit Greenland, clear-sky:
        %load sas_Smm_sfc_flx_frc_clr.txt;
        %sas_Smm_sfc_flx_frc_clr(find(sas_Smm_sfc_flx_frc_clr==0))=1e-30;
        %flx_slr = sas_Smm_sfc_flx_frc_clr;

        Fs(1:nbr_wvl,1) = flx_slr./(mu_not*pi);  % direct-beam incident flux
        Fd(1:nbr_wvl,1) = 0;                     % diffuse incident flux
    
    elseif (DIRECT == 0)
        % mid-latitude winter, cloudy-sky
        load mlw_sfc_flx_frc_cld.txt;
        mlw_sfc_flx_frc_cld(find(mlw_sfc_flx_frc_cld==0))=1e-30;
        flx_slr = mlw_sfc_flx_frc_cld;
        
        % Summit Greenland, cloudy:
        %load sas_Smm_sfc_flx_frc_cld.txt;
        %sas_Smm_sfc_flx_frc_cld(find(sas_Smm_sfc_flx_frc_cld==0))=1e-30;
        %flx_slr = sas_Smm_sfc_flx_frc_cld;

        Fd(1:nbr_wvl,1) = flx_slr; % direct-beam incident flux
        Fs(1:nbr_wvl,1) = 0;       % diffuse incident flux
  end;

  
  
elseif (BND_TYP == 2)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    %
    % SPECTRAL BANDS:
    %  1. 0.3-0.7 um
    %  2. 0.7-1.0 um
    %  3. 1.0-1.2 um
    %  4. 1.2-1.5 um
    %  5. 1.5-5.0 um

    % mid-latitude winter:
    if (DIRECT == 1)
        flx_slr(1,1) = 0.5028;
        flx_slr(2,1) = 0.2454;
        flx_slr(3,1) = 0.0900;
        flx_slr(4,1) = 0.0601;
        flx_slr(5,1) = 0.1017;
   
        Fs = flx_slr./(mu_not*pi);
        Fd(1:nbr_wvl,1) = 0;
  
    elseif (DIRECT == 0)
        flx_slr(1,1) = 0.5767;
        flx_slr(2,1) = 0.2480;
        flx_slr(3,1) = 0.0853;
        flx_slr(4,1) = 0.0462;
        flx_slr(5,1) = 0.0438;
        
        Fs(1:nbr_wvl,1) = 0;
        Fd = flx_slr;
    end;
    
    
elseif (BND_TYP == 3)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    %
    % SPECTRAL BANDS:
    %  1. 0.3-0.7 um
    %  2. 0.7-1.19 um
    %  3. 1.0-5.0 um
      
    if (DIRECT == 1)
        % MID-LATITUDE WINTER (SWNB2) (0.3-0.7, 0.7-1.19, 1.19-5.0)
        flx_slr(1,1) = 0.5028;
        flx_slr(2,1) = 0.3313;
        flx_slr(3,1) = 0.1659;
   
        Fs = flx_slr./(mu_not*pi);
        Fd(1:nbr_wvl,1) = 0;
  
    elseif (DIRECT == 0)
        flx_slr(1,1) = 0.5767;
        flx_slr(2,1) = 0.3297;
        flx_slr(3,1) = 0.0936;
        
        Fs(1:nbr_wvl,1) = 0;
        Fd = flx_slr;
    end;
    
end;



% Read Mie ice parameters from external NetCDF files (layer-specific):
for n=1:nbr_lyr
       
    %NEW code for reading in rds_coated() when rds_snw() is set to zero
       %(i.e. using ice spheres coated with a water film)
       
    if (rds_snw(n) == 0)    % JC EDIT
        
        %s3 = int2str(rds_coated(n));   % JC EDIT
        s3 = rds_coated(n);
        
        fl_in = strcat(wrkdir,fl_stb1,fl_stb3,s3,fl_stb2)   % JC EDIT
    %   *************************************************        
    elseif (rds_snw(n)<10)   % JC EDIT
        s1    = int2str(0);   % JC EDIT
        s2    = int2str(rds_snw(n));   % JC EDIT
        fl_in = strcat(wrkdir,fl_stb1,s1,s1,s1,s2,fl_stb2);   % JC EDIT
    
    elseif (rds_snw(n)<100)   % JC EDIT
        s1    = int2str(0);   % JC EDIT
        s2    = int2str(rds_snw(n));   % JC EDIT
        fl_in = strcat(wrkdir,fl_stb1,s1,s1,s2,fl_stb2);   % JC EDIT
    
    elseif ((rds_snw(n)>=100) & (rds_snw(n)<1000))   % JC EDIT
        s1    = int2str(0);   % JC EDIT
        s2    = int2str(rds_snw(n));   % JC EDIT
        fl_in = strcat(wrkdir,fl_stb1,s1,s2,fl_stb2);   % JC EDIT
    
    elseif ((rds_snw(n)>=1000) & (rds_snw(n)<9999))   % JC EDIT
        s1    = int2str(0);   % JC EDIT
        s2    = int2str(rds_snw(n));   % JC EDIT
        fl_in = strcat(wrkdir,fl_stb1,s2,fl_stb2);   % JC EDIT
    
    elseif (rds_snw(n)>9999)   % JC EDIT
        s2 = int2str(rds_snw(n));   % JC EDIT
        fl_in = strcat(wrkdir,fl_stb1,s2,fl_stb2);   % JC EDIT   
    
    elseif (rds_snw(n)==99999)   % JC EDIT
        fl_in = fl_in_usr;   % JC EDIT
    
    end;   % JC EDIT

    % single-scatter albedo, mass extinction coefficient, and
    % asymmetry paramater
    
    omega_snw(:,n)       = ncread(fl_in,'ss_alb');
    ext_cff_mss_snw(:,n) = ncread(fl_in,'ext_cff_mss');
    g_snw(:,n)           = ncread(fl_in,'asm_prm');
    
end


% Read Mie aerosol parameters (layer-independent)
fl_in1 = strcat(wrkdir,fl_sot1);
fl_in2 = strcat(wrkdir,fl_sot2);
fl_in3 = strcat(wrkdir,fl_dst1);
fl_in4 = strcat(wrkdir,fl_dst2);
fl_in5 = strcat(wrkdir,fl_dst3);
fl_in6 = strcat(wrkdir,fl_dst4);
fl_in7 = strcat(wrkdir,fl_ash1);
fl_in8 = strcat(wrkdir,fl_bio1); % JC EDIT
fl_in9 = strcat(wrkdir,fl_bio2); % JC EDIT
fl_in10 = strcat(wrkdir,fl_bio3); % JC EDIT
fl_in11 = strcat(wrkdir,fl_bio4); % JC EDIT
fl_in12 = strcat(wrkdir,fl_bio5); % JC EDIT
fl_in13 = strcat(wrkdir,fl_bio6); % JC EDIT
fl_in14 = strcat(wrkdir,fl_bio7); % JC EDIT
fl_in15 = strcat(wrkdir,fl_RBio1); % JC EDIT
fl_in16 = strcat(wrkdir,fl_hematite); % JC EDIT
fl_in17 = strcat(wrkdir,fl_mixed_sand); % JEC EDIT


omega_aer(:,1)       = ncread(fl_in1,'ss_alb');
g_aer(:,1)           = ncread(fl_in1,'asm_prm');
ext_cff_mss_aer(:,1) = ncread(fl_in1,'ext_cff_mss');
  
omega_aer(:,2)       = ncread(fl_in2,'ss_alb');
g_aer(:,2)           = ncread(fl_in2,'asm_prm');
ext_cff_mss_aer(:,2) = ncread(fl_in2,'ext_cff_mss_cor');
%ext_cff_mss_aer(:,2)= ncread(fl_in2,'ext_cff_mss');

omega_aer(:,3)       = ncread(fl_in3,'ss_alb');
g_aer(:,3)           = ncread(fl_in3,'asm_prm');
ext_cff_mss_aer(:,3) = ncread(fl_in3,'ext_cff_mss');

omega_aer(:,4)       = ncread(fl_in4,'ss_alb');
g_aer(:,4)           = ncread(fl_in4,'asm_prm');
ext_cff_mss_aer(:,4) = ncread(fl_in4,'ext_cff_mss');

omega_aer(:,5)       = ncread(fl_in5,'ss_alb');
g_aer(:,5)           = ncread(fl_in5,'asm_prm');
ext_cff_mss_aer(:,5) = ncread(fl_in5,'ext_cff_mss');

omega_aer(:,6)       = ncread(fl_in6,'ss_alb');   
g_aer(:,6)           = ncread(fl_in6,'asm_prm');  
ext_cff_mss_aer(:,6) = ncread(fl_in6,'ext_cff_mss');   

omega_aer(:,7)       = ncread(fl_in7,'ss_alb');   
g_aer(:,7)           = ncread(fl_in7,'asm_prm');  
ext_cff_mss_aer(:,7) = ncread(fl_in7,'ext_cff_mss'); 

omega_aer(:,8)       = ncread(fl_in8,'ss_alb');   % JC EDIT
g_aer(:,8)           = ncread(fl_in8,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,8) = ncread(fl_in8,'ext_cff_mss');   % JC EDIT

omega_aer(:,9)       = ncread(fl_in9,'ss_alb');   % JC EDIT
g_aer(:,9)           = ncread(fl_in9,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,9) = ncread(fl_in9,'ext_cff_mss');   % JC EDIT

omega_aer(:,10)       = ncread(fl_in10,'ss_alb');   % JC EDIT
g_aer(:,10)           = ncread(fl_in10,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,10) = ncread(fl_in10,'ext_cff_mss');   % JC EDIT

omega_aer(:,11)       = ncread(fl_in11,'ss_alb');   % JC EDIT
g_aer(:,11)           = ncread(fl_in11,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,11) = ncread(fl_in11,'ext_cff_mss');   % JC EDIT

omega_aer(:,12)       = ncread(fl_in12,'ss_alb');   % JC EDIT
g_aer(:,12)           = ncread(fl_in12,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,12) = ncread(fl_in12,'ext_cff_mss');   % JC EDIT

omega_aer(:,13)       = ncread(fl_in13,'ss_alb');   % JC EDIT
g_aer(:,13)           = ncread(fl_in13,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,13) = ncread(fl_in13,'ext_cff_mss');   % JC EDIT
 
omega_aer(:,14)       = ncread(fl_in14,'ss_alb');   % JC EDIT
g_aer(:,14)           = ncread(fl_in14,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,14) = ncread(fl_in14,'ext_cff_mss');    % JC EDIT

omega_aer(:,15)       = ncread(fl_in15,'ss_alb');   % JC EDIT
g_aer(:,15)           = ncread(fl_in15,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,15) = ncread(fl_in15,'ext_cff_mss');   % JC EDIT

omega_aer(:,16)       = ncread(fl_in16,'ss_alb');   % JC EDIT
g_aer(:,16)           = ncread(fl_in16,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,16) = ncread(fl_in16,'ext_cff_mss');   % JC EDIT

omega_aer(:,17)       = ncread(fl_in17,'ss_alb');   % JC EDIT
g_aer(:,17)           = ncread(fl_in17,'asm_prm');   % JC EDIT
ext_cff_mss_aer(:,17) = ncread(fl_in17,'ext_cff_mss');   % JC EDIT

% Set aerosol concentration matrix:
mss_cnc_aer(1:nbr_lyr,1) = mss_cnc_sot1;
mss_cnc_aer(1:nbr_lyr,2) = mss_cnc_sot2;
mss_cnc_aer(1:nbr_lyr,3) = mss_cnc_dst1;
mss_cnc_aer(1:nbr_lyr,4) = mss_cnc_dst2;
mss_cnc_aer(1:nbr_lyr,5) = mss_cnc_dst3;
mss_cnc_aer(1:nbr_lyr,6) = mss_cnc_dst4;
mss_cnc_aer(1:nbr_lyr,7) = mss_cnc_ash1;
mss_cnc_aer(1:nbr_lyr,8) = mss_cnc_bio1; % JC EDIT
mss_cnc_aer(1:nbr_lyr,9) = mss_cnc_bio2; % JC EDIT
mss_cnc_aer(1:nbr_lyr,10) = mss_cnc_bio3; % JC EDIT
mss_cnc_aer(1:nbr_lyr,11) = mss_cnc_bio4; % JC EDIT
mss_cnc_aer(1:nbr_lyr,12) = mss_cnc_bio5; % JC EDIT
mss_cnc_aer(1:nbr_lyr,13) = mss_cnc_bio6; % JC EDIT
mss_cnc_aer(1:nbr_lyr,14) = mss_cnc_bio7; % JC EDIT
mss_cnc_aer(1:nbr_lyr,15) = mss_cnc_RBio1; %JC EDIT
mss_cnc_aer(1:nbr_lyr,16) = mss_cnc_hematite; % JC EDIT
mss_cnc_aer(1:nbr_lyr,17) = mss_cnc_mixed_sand; % JC EDIT


% convert to units of kg/kg:
mss_cnc_aer = mss_cnc_aer.*10^-9;

% BEGIN RT SOLVER:

% Calculate effective tau, omega, g for the (snow+impurity) system
for n=1:nbr_lyr
    L_snw(n)     = rho_snw(n)*dz(n);                 % Snow column mass (kg/m^2) (array)
    tau_snw(:,n) = L_snw(n).*ext_cff_mss_snw(:,n);
    
    for j=1:nbr_aer
        L_aer(n,j)     = L_snw(n)*mss_cnc_aer(n,j);
        tau_aer(:,n,j) = L_aer(n,j).*ext_cff_mss_aer(:,j);
    end
  
    tau_sum(1:nbr_wvl,1)   = 0.0;
    omega_sum(1:nbr_wvl,1) = 0.0;
    g_sum(1:nbr_wvl,1)     = 0.0;
  
    for j=1:nbr_aer
        tau_sum   = tau_sum + tau_aer(:,n,j);
        omega_sum = omega_sum + (tau_aer(:,n,j).*omega_aer(:,j));
        g_sum     = g_sum + (tau_aer(:,n,j).*omega_aer(:,j).*g_aer(:,j));
    end
  
    tau(:,n)   = tau_sum + tau_snw(:,n);
    omega(:,n) = (1./tau(:,n)).*(omega_sum+ (omega_snw(:,n).*tau_snw(:,n)));
    g(:,n)     = (1./(tau(:,n).*omega(:,n))) .* (g_sum+ (g_snw(:,n).*omega_snw(:,n).*tau_snw(:,n)));
end
  
  
% Perform Delta-transformations, if called for
if (DELTA == 1)
    g_star     = g./(1+g);
    omega_star = ((1-(g.^2)).*omega) ./ (1-(omega.*(g.^2)));
    tau_star   = (1-(omega.*(g.^2))).*tau;
else
    g_star     = g;
    omega_star = omega;
    tau_star   = tau;
end;

% Calculate total column optical depth
% tau_clm(:,n) = total optical depth of layers above n, or optical
% depth from upper model boundary to top of layer n
tau_clm(1:nbr_wvl,1:1) = 0.0;
for n=2:nbr_lyr
    tau_clm(:,n) = tau_clm(:,n-1)+tau_star(:,n-1);
end


% Boundary Condition: Upward (reflected) direct beam flux at lower model boundary
% (Eq. 37)
S_sfc = R_sfc.*mu_not.*exp(-(tau_clm(:,nbr_lyr)+tau_star(:,nbr_lyr))./mu_not).*pi.*Fs;


% Apply 2-stream approximation technique (Toon et al., Table 1.)

if (APRX_TYP == 1)
    % Eddington:
    gamma1 = (7-(omega_star.*(4+(3*g_star))))./4;
    gamma2 = -(1-(omega_star.*(4-(3*g_star))))./4;
    gamma3 = (2-(3*g_star.*mu_not))./4;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
  
elseif (APRX_TYP == 2)
    % Quadrature: 
    gamma1 = sqrt(3).*(2-(omega_star.*(1+g_star)))./2;
    gamma2 = omega_star.*sqrt(3).*(1-g_star)./2;
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 1/sqrt(3);
  
elseif (APRX_TYP == 3)
    % Hemispheric mean:
    gamma1 = 2 - (omega_star.*(1+g_star));
    gamma2 = omega_star.*(1-g_star);
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
end;

% Eq. 21 and 22
lambda = sqrt(abs((gamma1.^2) - (gamma2.^2)));
GAMMA  = gamma2./(gamma1+lambda);

% Eq. 44
e1 = 1+(GAMMA.*exp(-lambda.*tau_star));
e2 = 1-(GAMMA.*exp(-lambda.*tau_star));
e3 = GAMMA + exp(-lambda.*tau_star);
e4 = GAMMA - exp(-lambda.*tau_star);


% Before calculating the C functions:
% C functions are indeterminate if [lambda^2 = 1/(mu_not^2)]
% upper bound of lambda^2 = 4 for delta-hemispheric mean, and
% lambda^2 = 3 for delta-Eddington.  So, problem can only arise
% when 0.5 < mu_not < 1.0.  
% Assuming that abs(lambda^2 - 1/(mu_not^2)) < 0.01 is dangerous:
% Let x= lambda^2 - 1/(mu_not^2),
% dx/dmu_not = 2*mu_not^-3
% dx = 0.01 = 2*(mu_not^-3)*dmu_not
% For range of possible mu_not:
% 0.04 > dmu_not > 0.005
% So, changing mu_not by 0.04 will be sufficient:

% Note: not implemented here because of vectorized code, 
% but an example is shown.  MUST be implemented in CLM

% Calculate C function for each layer
% (evaluated at the top and bottom of each layer n)
% (Eq. 23 and 24)

for n=1:nbr_lyr

    if (sum(Fs) > 0.0)
      
        % Stability check, see explanation above:
        %for p=1:nbr_wvl
        %temp1 = ( lambda(p,n)^2) - 1/mu_not^2);
        %if abs(temp) < 0.01
        %  mu_not_sv = mu_not;
        %  if (temp1 > 0)
        %    mu_not = mu_not + 0.04;
        %  else
        %    mu_not = mu_not - 0.04
        %  end;
        %end;
        %end
        %if (n==1)
        %    mu_not=mu_not+0.00;
        %else
        %    mu_not = mu_not-0.00;
        %end;
	
        C_pls_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+ ...
                           (gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
  
            
        C_mns_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+ ...
                           (gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));


        C_pls_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+(gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
            
        C_mns_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+(gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
    else
        % no direct-beam flux:
        C_pls_btm(1:nbr_wvl,n) = 0.0;
        C_mns_btm(1:nbr_wvl,n) = 0.0;
        C_pls_top(1:nbr_wvl,n) = 0.0;
        C_mns_top(1:nbr_wvl,n) = 0.0;
    end;
      
end
  

% Eq. 41-43
for i=1:2*nbr_lyr
    % Boundary values for i=1 and i=2nbr_lyr, specifics for i=odd and i=even    
    if (i==1)
        A(1:nbr_wvl,1) = 0.0;
        B(1:nbr_wvl,1) = e1(:,1);
        D(1:nbr_wvl,1) = -e2(:,1);
        E(1:nbr_wvl,1) = Fd(:)-C_mns_top(:,1);
					  
    elseif(i==2*nbr_lyr)
        A(:,i) = e1(:,nbr_lyr)-(R_sfc.*e3(:,nbr_lyr));
        B(:,i) = e2(:,nbr_lyr)-(R_sfc.*e4(:,nbr_lyr));
        D(:,i) = 0.0;
        E(:,i) = S_sfc(:) - C_pls_btm(:,nbr_lyr) + (R_sfc.*C_mns_btm(:,nbr_lyr));
    
    elseif(mod(i,2)==1)  % If odd and i>=3 (n=1 for i=3. this is confusing)
        n      = floor(i/2);
        A(:,i) = (e2(:,n).*e3(:,n))-(e4(:,n).*e1(:,n));
        B(:,i) = (e1(:,n).*e1(:,n+1))-(e3(:,n).*e3(:,n+1));
        D(:,i) = (e3(:,n).*e4(:,n+1))-(e1(:,n).*e2(:,n+1));
        E(:,i) = (e3(:,n).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e1(:,n).*(C_mns_btm(:,n)-C_mns_top(:,n+1)));
    
    elseif(mod(i,2)==0)  % If even and i<=2nbr_lyr
        n      = (i/2);
        A(:,i) = (e2(:,n+1).*e1(:,n))-(e3(:,n).*e4(:,n+1));
        B(:,i) = (e2(:,n).*e2(:,n+1))-(e4(:,n).*e4(:,n+1));
        D(:,i) = (e1(:,n+1).*e4(:,n+1))-(e2(:,n+1).*e3(:,n+1));
        E(:,i) = (e2(:,n+1).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e4(:,n+1).*(C_mns_top(:,n+1)-C_mns_btm(:,n))); 
    end;
end


% Eq. 45
AS(1:nbr_wvl,2*nbr_lyr) = A(:,2*nbr_lyr)./B(:,2*nbr_lyr);
DS(1:nbr_wvl,2*nbr_lyr) = E(:,2*nbr_lyr)./B(:,2*nbr_lyr);
  
% Eq. 46
for i=(2*nbr_lyr-1):-1:1
    X(1:nbr_wvl,i) = 1./(B(:,i)-(D(:,i).*AS(:,i+1)));
    AS(:,i)        = A(:,i).*X(:,i);
    DS(:,i)        = (E(:,i)-(D(:,i).*DS(:,i+1))).*X(:,i);
end

% Eq. 47
Y(1:nbr_wvl,1) = DS(:,1);
for i=2:2*nbr_lyr
    Y(:,i) = DS(:,i)-(AS(:,i).*Y(:,i-1));
end;


for n=1:nbr_lyr
    % Direct beam flux at the base of each layer (Eq. 50)
    direct(1:nbr_wvl,n) = mu_not*pi*Fs.*exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not);

    % Net flux (positive upward = F_up-F_down) at the base of each
    % layer (Eq. 48)
    F_net(1:nbr_wvl,n) = (Y(:,(2*n-1)).*(e1(:,n)-e3(:,n))) +...
        (Y(:,(2*n)).*(e2(:,n)-e4(:,n))) + ...
        C_pls_btm(:,n) - C_mns_btm(:,n) - direct(:,n);

    % Mean intensity at the base of each layer (Eq. 49):
    intensity(1:nbr_wvl,n) = (1/mu_one).*...
        ( Y(:,(2*n-1)).*(e1(:,n)+e3(:,n)) + ...
          Y(:,(2*n)).*(e2(:,n)+e4(:,n)) + C_pls_btm(:,n) + C_mns_btm(:,n)) +...
        (direct(:,n)./mu_not);
  
    intensity(1:nbr_wvl,n) = intensity(1:nbr_wvl,n)./(4*pi);
end
  

% Upward flux at upper model boundary (Eq. 31):
F_top_pls = (Y(:,1).*(exp(-lambda(:,1).*tau_star(:,1))+GAMMA(:,1))) + ...
    (Y(:,2).*(exp(-lambda(:,1).*tau_star(:,1))-GAMMA(:,1))) + ...
    C_pls_top(:,1);


for n=1:nbr_lyr
    % Upward flux at the bottom of each layer interface (Eq. 31)
    F_up(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(exp(0) + GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        Y(:,2*n)  .*(exp(0) - GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        C_pls_btm(:,n);
  
    % Downward flux at the bottom of each layer interface (Eq. 32,
    % plus direct-beam component):
    F_down(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(GAMMA(:,n).*exp(0) + exp(-lambda(:,n).*tau_star(:,n))) + ...
        Y(:,2*n)  .*(GAMMA(:,n).*exp(0) - exp(-lambda(:,n).*tau_star(:,n))) + ...
        C_mns_btm(:,n) + direct(:,n);

    % Derived net (upward-downward) flux
    % (should equal F_net, Eq. 48)
    F_net2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) - F_down(1:nbr_wvl,n);
  
    % planar intensity:
    intensity2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) + F_down(1:nbr_wvl,n);
end

% surface planar intensity:
intensity2_top(1:nbr_wvl) = F_top_pls + ((mu_not*pi*Fs)+Fd);

% diagnostic:
if (BND_TYP==1)
    intensity_out(1)           = intensity2_top(26);
    intensity_out(2:nbr_lyr+1) = intensity2(26,1:nbr_lyr);
end;


% Net flux at lower model boundary = bulk transmission through entire
% media = absorbed radiation by underlying surface:
F_btm_net = -F_net(:,nbr_lyr);


% Hemispheric wavelength-dependent albedo:
if (BND_TYP < 4)
    albedo = F_top_pls./((mu_not*pi*Fs)+Fd);
end


% Net flux at upper model boundary
F_top_net(1:nbr_wvl,1) = F_top_pls - ((mu_not*pi*Fs)+Fd);

  
% Absorbed flux in each layer (negative if there is net emission (bnd_typ==4))
for n=1:nbr_lyr
    if(n==1)
        F_abs(1:nbr_wvl,1) = F_net(:,1)-F_top_net;
    else
        F_abs(:,n) = F_net(:,n) - F_net(:,n-1);
    end;
end


% Set indices for VIS and NIR
if (BND_TYP == 1)
    vis_max_idx = 40;
    nir_max_idx = length(wvl);
elseif ((BND_TYP == 2) | (BND_TYP == 3))
    vis_max_idx = 1;
    nir_max_idx = length(wvl);
end;


% Spectrally-integrated absorption in each layer:
F_abs_slr=sum(F_abs);
for n=1:nbr_lyr
    F_abs_vis(n) = sum(F_abs(1:vis_max_idx,n));
    F_abs_nir(n) = sum(F_abs(vis_max_idx+1:nir_max_idx,n));
end

% Spectrally-integrated absorption by underlying surface:
F_abs_btm = sum(F_btm_net);
F_abs_vis_btm = sum(F_btm_net(1:vis_max_idx));
F_abs_nir_btm = sum(F_btm_net(vis_max_idx+1:nir_max_idx));


% Radiative heating rate:
heat_rt = F_abs_slr./(L_snw.*2117);   % [K/s] 2117 = specific heat ice (J kg-1 K-1)	   
heat_rt = heat_rt.*3600;              % [K/hr]
				      
				      
% Energy conservation check:
% Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
energy_sum = (mu_not*pi*Fs)+Fd - (sum(F_abs,2) + F_btm_net + F_top_pls);

if (sum(abs(energy_sum)) > 1e-10)
    energy_conservation_error = sum(abs(energy_sum))
    %error(strcat('Energy conservation error of: ',num2str(sum(abs(energy_sum)))));
end;

% spectrally-integrated terms (remove semi-colons to write-out
% these values):
sum(energy_sum);        % energy conservation total error
sum((mu_not*pi*Fs)+Fd); % total incident insolation (W m-2)
sum(sum(F_abs));        % total energy absorbed by all snow layers
sum(F_btm_net);         % total energy absorbed by underlying substrate


% Spectrally-integrated solar, visible, and NIR albedos:
alb     = sum(flx_slr.*albedo)./sum(flx_slr);

alb_vis = sum(flx_slr(1:vis_max_idx).*albedo(1:vis_max_idx))/...
          sum(flx_slr(1:vis_max_idx));

alb_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*albedo(vis_max_idx+1:nir_max_idx))/...
          sum(flx_slr(vis_max_idx+1:nir_max_idx));


% Diagnostic for comparing 470-band solutions with 5-band solutions:
%  Spectrally-integrated 5-band albedo:
if (BND_TYP == 1)
    bnd1a=1;
    bnd1b=40;
    bnd2a=41;
    bnd2b=70;
    bnd3a=71;
    bnd3b=90;
    bnd4a=91;
    bnd4b=120;
    bnd5a=121;
    bnd5b=470;
  
    bnd6a=1;
    bnd6b=40;
    bnd7a=41;
    bnd7b=89;
    bnd8a=90;
    bnd8b=470;
  
    alb1 = sum(flx_slr(bnd1a:bnd1b).*albedo(bnd1a:bnd1b))/sum(flx_slr(bnd1a:bnd1b));
    alb2 = sum(flx_slr(bnd2a:bnd2b).*albedo(bnd2a:bnd2b))/sum(flx_slr(bnd2a:bnd2b));
    alb3 = sum(flx_slr(bnd3a:bnd3b).*albedo(bnd3a:bnd3b))/sum(flx_slr(bnd3a:bnd3b));
    alb4 = sum(flx_slr(bnd4a:bnd4b).*albedo(bnd4a:bnd4b))/sum(flx_slr(bnd4a:bnd4b));
    alb5 = sum(flx_slr(bnd5a:bnd5b).*albedo(bnd5a:bnd5b))/sum(flx_slr(bnd5a:bnd5b));
    alb6 = sum(flx_slr(bnd6a:bnd6b).*albedo(bnd6a:bnd6b))/sum(flx_slr(bnd6a:bnd6b));
    alb7 = sum(flx_slr(bnd7a:bnd7b).*albedo(bnd7a:bnd7b))/sum(flx_slr(bnd7a:bnd7b));
    alb8 = sum(flx_slr(bnd8a:bnd8b).*albedo(bnd8a:bnd8b))/sum(flx_slr(bnd8a:bnd8b));
end;

% Spectrally-integrated VIS and NIR total snowpack absorption:
abs_vis = sum(flx_slr(1:vis_max_idx).*(1-albedo(1:vis_max_idx)));
abs_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*(1-albedo(vis_max_idx+1:nir_max_idx)));



% Output diagnostics:
%  1. The basics:
data_out(:,1) = wvl;            % spectral wavelength bands (um)
data_out(:,2) = albedo;         % spectral hemispheric albedo
data_out(1,3) = alb;            % solar broadband albedo
data_out(2,3) = alb_vis;        % visible (0.3-0.7um) albedo
data_out(3,3) = alb_nir;        % near-IR (0.7-5.0um) albedo
data_out(4,3) = sum(F_abs_slr); % total radiative absorption by all snow layers (not including underlying substrate)
data_out(5,3) = F_abs_slr(1);   % top layer solar absorption
if (nbr_lyr > 1)
    data_out(6,3) = F_abs_slr(2); % 2nd layer solar absorption
else
    data_out(6,3) = NaN;          % 2nd layer solar absorption
end;
data_out(:,5)  = sum(F_abs,2);  % spectral absorption

% more detail:
if (BND_TYP == 1)
    % different band-weighted albedos:
    data_out(1,4) = alb1;
    data_out(2,4) = alb2;
    data_out(3,4) = alb3;
    data_out(4,4) = alb4;
    data_out(5,4) = alb5;
    data_out(7,3) = alb6;
    data_out(8,3) = alb7;
    data_out(9,3) = alb8;
end;

data_out(6,4) = F_abs_slr(1); % top layer solar absorption
data_out(7,4) = F_abs_vis(1); % top layer VIS absorption
data_out(8,4) = F_abs_nir(1); % top layer NIR absorption

if (nbr_lyr > 1)
    data_out(9,4) = F_abs_slr(2);
    data_out(10,4) = F_abs_vis(2);
    data_out(11,4) = F_abs_nir(2);
    if (nbr_lyr > 2)
        data_out(12,4) = F_abs_slr(3);
        data_out(13,4) = F_abs_vis(3);
        data_out(14,4) = F_abs_nir(3);
        if (nbr_lyr > 3)
            data_out(15,4) = F_abs_slr(4);
            data_out(16,4) = F_abs_vis(4);
            data_out(17,4) = F_abs_nir(4);
            if (nbr_lyr > 4)
                data_out(18,4) = F_abs_slr(5);
                data_out(19,4) = F_abs_vis(5);
                data_out(20,4) = F_abs_nir(5);
            end;
        end;
    end;
end;

data_out(18,4) = F_abs_btm;     % solar absorption by underlying surface
data_out(19,4) = F_abs_vis_btm; % VIS absorption by underlying surface
data_out(20,4) = F_abs_nir_btm; % NIR absorption by underlying surface
data_out(21,4) = sum((mu_not*pi*Fs))+sum(Fd);  % total downwelling energy on upper boundary
data_out([1:5],6) = heat_rt; % JC EDIT output radiative heating rate in each layer in K/hr 
data_out([1:470],[7:11]) = intensity2; %JC EDIT output planar intensity per layer




% plot modeled albedo:
if (1==0)
    plot(wvl,albedo,'k-');
    axis([0.3 2.5 0 1]);
    grid on;
end;

