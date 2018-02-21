
function [alb_slr,albedo_spectral] = SNICAR_function_CASPA(BND_TYP,DIRECT,APRX_TYP,DELTA,coszen,R_sfc,dz,rho_snw,rds_snw,x)

nbr_lyr  = length(dz);  % number of snow layers
nbr_aer = 16;


% PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)
mss_cnc_sot1(1:nbr_lyr)  = [0,0,0,0,0];  % uncoated black carbon
mss_cnc_sot2(1:nbr_lyr)  = [0,0,0,0,0];    % coated black carbon
mss_cnc_dst1(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 1
mss_cnc_dst2(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 2
mss_cnc_dst3(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 3
mss_cnc_dst4(1:nbr_lyr)  = [0,0,0,0,0];    % dust species 4
mss_cnc_ash1(1:nbr_lyr)  = [0,0,0,0,0];    % volcanic ash species 1
mss_cnc_bio1(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 1
mss_cnc_bio2(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 2
mss_cnc_bio3(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 3
mss_cnc_bio4(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 4
mss_cnc_bio5(1:nbr_lyr)  = [x,0,0,0,0];    % Biological impurity species 5
mss_cnc_bio6(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 6
mss_cnc_bio7(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 6
mss_cnc_water1(1:nbr_lyr) = [0,0,0,0,0];   % Water, 2 mm spheres
mss_cnc_hematite(1:nbr_lyr) = [0,0,0,0,0];   % Water, 2 mm spheres
mss_cnc_mixed_sand(1:nbr_lyr) = [0,0,0,0,0];


% FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
fl_dst1  = 'aer_dst_bln_20060904_01.nc';
fl_dst2  = 'aer_dst_bln_20060904_02.nc';
fl_dst3  = 'aer_dst_bln_20060904_03.nc';
fl_dst4  = 'aer_dst_bln_20060904_04.nc';
fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
fl_bio1  = 'biological_1.nc'; % Biological impurity 1 (30um diameter, 1.5%chll a,10% each 1 & 2 carotenoids) )
fl_bio2  = 'biological_2.nc'; % Biological impurity 2 (30um diameter, 1.5%chll a, 5% each 1 % 2 carotenoids)
fl_bio3  = 'biological_3.nc'; % Biological impurity 3 (30um diameter, 1.5%chll a, 1% each 1 % 2 carotenoids)
fl_bio4  = 'biological_4.nc'; % Biological impurity 4 (30um diameter, 1.5%chll a only)
fl_bio5  = 'biological_5.nc'; % Biological impurity 5 (10um diameter, pigs as per bio2)
fl_bio6  = 'biological_6.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
fl_bio7 = 'biological_7.nc';
fl_hematite  = 'Hematite.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
fl_mixed_sand  = 'mixed_sand.nc'; % Mixed sand (quartz and clays)
rds_coated = [0,0,0,0,0];

% call SNICAR with these inputs:
data_in = snicar8d(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
    dz, rho_snw, rds_snw, rds_coated, nbr_aer, mss_cnc_sot1, ...
    mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
    mss_cnc_dst3, mss_cnc_dst4, mss_cnc_ash1, mss_cnc_bio1, mss_cnc_bio2,mss_cnc_bio3,mss_cnc_bio4,mss_cnc_bio5, mss_cnc_bio6,mss_cnc_bio7, mss_cnc_hematite, mss_cnc_mixed_sand, fl_sot1, ...
    fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4, fl_ash1, fl_bio1,fl_bio2,fl_bio3,fl_bio4,fl_bio5,fl_bio6,fl_bio7, fl_hematite, fl_mixed_sand);


% process input data:
% (see description of data_out at the end of snicar8d.m for more info)
wvl         = data_in(:,1);   % wavelength grid
albedo      = data_in(:,2);   % spectral albedo
alb_slr     = data_in(1,3);   % broadband albedo (0.3-5.0um)
alb_vis     = data_in(2,3);   % visible albedo (0.3-0.7um)
alb_nir     = data_in(3,3);   % near-IR albedo (0.7-5.0um)
flx_abs_snw = data_in(4,3);   % total radiative absorption by all snow layers (not including underlying substrate)

flx_abs(1)     = data_in(6,4); % top layer solar absorption
flx_vis_abs(1) = data_in(7,4); % top layer VIS absorption
flx_nir_abs(1) = data_in(8,4); % top layer NIR absorption

alb_slr;
albedo_spectral = albedo;

  
end
    
