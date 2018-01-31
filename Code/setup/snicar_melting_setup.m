% Driver routine for SNICAR.  Also see commenting in snicar8d.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% JC EDIT: This version is designed to interface with the CASPA cellular
% automaton for predicting albedo change over space and time. This version
% predicts the grain growth resulting from heating of each layer of the
% snowpack. Running this driver outputs the spectral and broadband
% albedo and also the new grain sizes after 1 CASPA timestep. 

% This is run from the driver script snicar_melting_driver.m, which also
% calls the mie solvers to ensure the right mie optical property files are
% available in the working directory for each grain produced by this code.

% This code accounts for wet and dry snow ageing according to the equations
% presented by Flanner and Zender (2006) and implemented in CLM. Some
% additions are also made to account for water. Water is generated in each
% layer according to the ecess energy above freezing point. All this energy
% is assumed to result in melting a certain mass of water. The mass of ice
% in the layer is reduced accordingly. The ratio of ice to water mas sis
% used to calculate the water fraction. This is then used to determine wet
% grain ageing in the layer. The water can then percolate downwards when it
% has accumulated past a certain threshold, provided the lower layer is not
% saturated. If liquid water arrives in a layer and the layer temperature
% is sub-freezing, it becomes refrozen ice with Reff = 1500. The new Reff
% of the layer is then the average of the aged grain radius and the
% refrozen ice radius weighted by the refrozen ice fraction. Where thr
% temperature is at or above freezing, the liquid water mass is converted
% into a thickness of water coating around the grain. Grains of appropriate
% size and with a specific water coating can then be modelled using the mie
% solver



function [overall_new_r, f_refs, new_T, fliqs, dz] = snicar_melting_setup(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, rds_snw, rds_coated, x, initial_T, fliqs, f_refs)

nbr_lyr  = length(dz);  % number of snow layers

% IF COATED GRAINS USED, SET rds_snw() to ZEROS and use rds_coated()
% IF UNCOATED GRAINS USED, SET rds_coated to ZEROS and use rds_snw()


% NUMBER OF AEROSOL SPECIES IN SNOW (ICE EXCLUDED)
%  Species numbers (used in snicar8d.m) are:
%    1: uncoated black carbon
%    2: coated black carbon
%    3: dust size 1
%    4: dust size 2
%    5: dust size 3
%    6: dust size 4
%    7: volcanic ash
%    8: biological impurity 1
%    9: biological impurity 2
%    10: biological impurity 3
%    11: biological impurity 4
%    12: biological impurity 5
%    13: biological impurity 6
%    14: biological impurity 7
%    15: hematite
%    16: mixed sand (quartz and clay)

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
mss_cnc_bio2(1:nbr_lyr)  = [x,0,0,0,0];    % Biological impurity species 2
mss_cnc_bio3(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 3
mss_cnc_bio4(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 4
mss_cnc_bio5(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 5
mss_cnc_bio6(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 6
mss_cnc_bio7(1:nbr_lyr)  = [0,0,0,0,0];    % Biological impurity species 7
mss_cnc_hematite(1:nbr_lyr) = [0,0,0,0,0];   % hematite
mss_cnc_mixed_sand(1:nbr_lyr) = [0,0,0,0,0]; % mixed quartz/clay sand


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
fl_bio7  = 'biological_7.nc'; % Biological impurity 7 (20um diameter, pigs as per bio2)
fl_hematite  = 'Hematite.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
fl_mixed_sand  = 'mixed_sand.nc'; % Mixed sand (quartz and clays)

% call SNICAR with these inputs:
data_in = snicar8d(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
    dz, rho_snw, rds_snw, rds_coated, nbr_aer, mss_cnc_sot1, ...
    mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
    mss_cnc_dst3, mss_cnc_dst4, mss_cnc_ash1, mss_cnc_bio1, mss_cnc_bio2,mss_cnc_bio3,mss_cnc_bio4,mss_cnc_bio5, mss_cnc_bio6, mss_cnc_bio7, mss_cnc_hematite, mss_cnc_mixed_sand, fl_sot1, ...
    fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4, fl_ash1, fl_bio1,fl_bio2,fl_bio3,fl_bio4,fl_bio5,fl_bio6, fl_bio7, fl_hematite, fl_mixed_sand);


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
heat_rt = data_in([1:5],6); % heating rate per layer K/hr
F_abs = [data_in(6,4),data_in(9,4),data_in(12,4),data_in(15,4),data_in(18,4)]; % absorbed energy per layer (W/m2)
%albedo = smooth(albedo,0.005); % add a simple smoothing function with short period

%Report albedo

alb_slr; % albedo over solar spectrum
flx_abs_snw; % absorbed energy in the snowpack




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SNOW GRAIN EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section of code uses some SNICAR outputs, some user defined
% variables and a call to a function (grain_size_evolution) in a separate
% script (grain_evolution.m) saved in the working directory. The separate
% script requires access to a lookup table (drdt_bst_fit_100.mat) saved in
% the working directory. This snow grain evolution calculation is an
% implementation of the parameterisation described by Flanner and Zender
% (2006) and used in the CLM model
% (http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% VARIABLE ASSIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User Defined

perc_threshold = 0.01; % max liquid water fraction that can be held by a layer without initiating percolation
doubling_time = 3; % must be the same as in CASPA (default = 3)
r0 = 100; % reff of fresh snow
i_mass = dz .* rho_snw; % (m * kgm-1)
w_mass = dz .* rho_snw.*fliqs; % (m * kgm-1)
underlying_temp = 260; % temp of underlying ice/ground

new_r = [0,0,0,0,0]; % set up empty array for new grain radius
new_T = [0,0,0,0,0]; % set up empty array for new temperature profile


% VARIABLES FROM SNICAR

heat_rt; % radiative heating rate in K/hr
heat_day = heat_rt.*24; % heating over day
heat_timestep = heat_day * 3; % heating over timestep (1 timestep = 3 days)
heat_timestep = heat_timestep.';
new_T = initial_T + heat_timestep;

% CONSTANTS

CICE = 2.11727e3; % Specific heat capacity of ice J/kg/K
CWAT = 4.188e3; % Specific heat capacity of water J/kg/K
Lf = 3.337e5; % Latent heat of fusion


%%%%%%%%%%%%%%%%%%%%%%% CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:nbr_lyr   % iterate through each vertical layer
    
    if i <= nbr_lyr-1
     
        temp_grad(i) = (new_T(i) - new_T(i+1)) / dz(i); % temp grad between successive snow layers
    else
        
        temp_grad(i) = (new_T(i) - underlying_temp) / dz(i); % temp grad between bottom layer and underlying ice/ground
    end   
    
    new_r(i) = grain_size_evolution(new_T(i),temp_grad(i),rho_snw(i),r0,rds_snw(i),doubling_time,fliqs(i),f_refs(i)); % function call from grain_size_evolution.m
    
        if new_T(i) >= 273.15  % if layer temp is greater than freezing point, initiate melting code to update liquid water fraction per layer...
            
            T_excess(i) = new_T(i) - 273.15; % calculate residual T that could be used for melting          
            c(i) = (w_mass(i) / dz(i) * CWAT) +  (i_mass(i) / dz(i) * CICE); % c is the volumetric ice/water heat capacity (J m-3 K-1)
            E_excess(i) = T_excess(i) * c(i); % excess energy (J/m3) = excess heat (K)* volumetric specific heat capacity of ice/water mix (J/m3/K) (OVER TIMESTEP)
            W_change(i) = (E_excess(i) / Lf)*dz(i); % change in ice mass = excess energy / energy required to melt ice (latent heat of fusion) * layer thickness (proxy for vol ice)
            new_w_mass(i) = w_mass(i) + W_change(i); % updated water mass = old + change
            i_mass(i) = i_mass(i) - W_change(i); % ice converted to water, so subtract change in water mass from initial ice mass
            fliqs(i) = new_w_mass(i) / i_mass(i); % new liquid water fraction

            % using this scheme assumes all excess energy is dissipated in
            % phase change ice -> water. Amount of excess energy determines
            % mass of water produced
            
        end              
end


% The following loop percolates new liquid water down through the snowpack.
% Regardless of the layer the water forms in, it will immediately percolate
% to the bottom, then accumulate from bottom up. There is a percolation
% threshold. If the liquid water fraction exceeds that threshold, liquid
% water percolates downwards, provided the layer beneath its max capacity 
% (fliq_max). 

for i = 1:1:nbr_lyr
    % if higher layer contains more water than lower layer, higher layer
    % exceeds percolation threshold and lower layer is not completely
    % saturated with liquid of refrozen water, percolate everything over
    % the percolation threshold down one layer (can cascade to btm layer)
    
    if i < nbr_lyr
        
        if fliqs(i) > fliqs(i+1) && fliqs(i) > perc_threshold && f_refs(i+1) < 1 
        
            fliq_diff = fliqs(i) - perc_threshold;
            fliqs(i+1) = fliqs(i+1) + fliq_diff;
            fliqs(i) = perc_threshold; 
        
        else
        
            fliqs(i) = fliqs(i); % if criteria not satisfied, water remains in layer
        end
        
    end % set bottom layer to max 1 liquid water
    
    % make sure liquid water fraction is between 0 and 1
    
    if fliqs(i) > 1
        
        fliqs(i) = 1;
    
    end
    
    if fliqs(i) < 0
        
        fliqs(i) = 0; 
    end
  
    
% Check whether liquid water exists in layers where temp is below freezing.
% If so, turn this liquid water fraction into refrozen fraction.    
   
    if fliqs(i) > 0 && new_T(i) < 273.15   % if there is water and the temp < freezing
        
        f_refs(i) = f_refs(i) + fliqs(i); % add liquid water fraction to refrozen fraction    
        
        fliqs(i) = 0; % remove liquid water from layer (it has moved to refrozen layer)
        
    end
        
    if f_refs(i) > 1
        
        f_refs(i) = 1;
    
    end 
    
    % Update layer thickness according to mass
        
    dz(i) = i_mass(i) / rho_snw(i); % update layer thickness (kg / kg m-1 = m)
    
    % Ensure layer does not grow vertically when it accrues mass - update
    % density instead
    
    if dz(i) > 0.05    
        
        rho_snw(i) = i_mass(i)/dz(i);
        dz(i) = 0.05;
    
    end    
    
end

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CHECKS AND BALANCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the layer must be described using one representative effective radius.
% Since we assume refrozen water to have reff = 1500um, the overall radius
% is the weighted average of the new radius and the refrozen grains
% (weighted by new_f_refs)

weight_r = 1-f_refs;

overall_new_r = (new_r.*weight_r) + (f_refs.*1500);

% since refrozen ice is assumed to have Reff = 1500, make this max Reff
% possible in any layer

for i = 1:1:nbr_lyr
    
    if (overall_new_r(i) > 1500) && (overall_new_r(i) <= 1750)
        overall_new_r(i) = 1500;
    end

% do not allow layer to be thinner then grain diameter

    if dz(i) < 2* overall_new_r(i)*10e-6
        FLAG = "R > LAYER"
        dz(i) = 2 * overall_new_r(i) * 10e-6; % layer thickness can't be less than grain diameter (2*r).
    end
    
% Do not allow temperature to exceed freezing point - this energy has
% already been dissipated by melting

    if new_T(i) > 273.15
        new_T(i) = 273.15; % reset max temps to freezing point
    end
    
end



end


