% Driver routine for SNICAR.  Also see commenting in snicar8d.m

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
% mss_cnc_bio1: mass mixing ratio of biological impurity species 1
%                  (units of parts per billion, ng g-1)
% mss_cnc_water1: mass mixing ratio of water type 1
%                  (units of parts per billion, ng g-1)

% fl_sot1:      name of file containing optical properties for BC species 1
% fl_sot2:      name of file containing optical properties for BC species 2
% fl_dst1:      name of file containing optical properties for dust species 1
% fl_dst2:      name of file containing optical properties for dust species 2
% fl_dst3:      name of file containing optical properties for dust species 3
% fl_dst4:      name of file containing optical properties for dust species 4
% fl_ash1:      name of file containing optical properties for ash species 1
% fl_bio1:      name of file containing optical properties for bio species 1
% fl_water1:    name of file containing optical properties for water type 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% RADIATIVE TRANSFER CONFIGURATION:
BND_TYP  = 1;        % 1= 470 spectral bands
DIRECT   = 1;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 3;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta

% COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM
coszen   = 0.50;

% REFLECTANCE OF SURFACE UNDERLYING SNOW:
%   Value is applied to all wavelengths.
%   User can also specify spectrally-dependent ground albedo
%   internally in snicar8d.m
R_sfc    = 0.15;


% SNOW LAYER THICKNESSES (array) (units: meters):
dz       = [0.01 0.02 0.02 0.02 0.02];
 
nbr_lyr  = length(dz);  % number of snow layers

% SNOW DENSITY OF EACH LAYER (units: kg/m3)
rho_snw(1:nbr_lyr) = [150, 150, 150, 150, 150];  


% SNOW EFFECTIVE GRAIN SIZE FOR EACH LAYER (units: microns):

rds_snw(1:nbr_lyr) = [1000, 1000, 1500, 2000, 2000];

% IF COATED GRAINS USED, SET rds_snw() to ZEROS and use rds_coated()
% IF UNCOATED GRAINS USED, SET rds_coated to ZEROS and use rds_snw()
% rds_coated should be input as strings.
rds_coated(1:nbr_lyr) = ["0","0","0","0","0"];


nbr_aer = 17;


for x = [0, 1e5 1e6 2e6 5e6 8e6 10e6]   % for reference: 1e3 = 1ug/g (1000 ppb or 1 ppm)
                  % 1 e6 = 1000ug = 1mg

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
mss_cnc_RBio1(1:nbr_lyr) = [0,0,0,0,0]; % Realistic Cell (measured pigments, 20 micron diameter)
mss_cnc_hematite(1:nbr_lyr) = [0,0,0,0,0];   % hematite
mss_cnc_mixed_sand(1:nbr_lyr) = [0,0,0,0,0];  % mixed sand


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
fl_RBio1 = 'Real_Bio_1.nc'; % Biological impurity with measure pigments, 20 micron diam
fl_hematite  = 'Hematite.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
fl_mixed_sand  = 'mixed_sand.nc'; % Mixed sand (quartz and clays)

% call SNICAR with these inputs:

data_in = snicar8d(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
    dz, rho_snw, rds_snw, rds_coated, nbr_aer, mss_cnc_sot1, ...
    mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
    mss_cnc_dst3, mss_cnc_dst4, mss_cnc_ash1, mss_cnc_bio1, mss_cnc_bio2,mss_cnc_bio3,mss_cnc_bio4,mss_cnc_bio5, mss_cnc_bio6, mss_cnc_bio7, mss_cnc_RBio1,mss_cnc_hematite, mss_cnc_mixed_sand, fl_sot1, ...
    fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4, fl_ash1, fl_bio1,fl_bio2,fl_bio3,fl_bio4,fl_bio5,fl_bio6, fl_bio7, fl_RBio1, fl_hematite, fl_mixed_sand);


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
intensity2 = data_in([1:470],[7:11]);
%albedo = smooth(albedo,0.005); % add a simple smoothing function with short period

%separate subsurface light field into layers
sub1 = smooth(intensity2(:,1),4);
sub2 = smooth(intensity2(:,2),4);
sub3 = smooth(intensity2(:,3),4);
sub4 = smooth(intensity2(:,4),4);
sub5 = smooth(intensity2(:,5),4);

sub1_tot = sum(sub1);
sub2_tot = sum(sub2);
sub3_tot = sum(sub3);
sub4_tot = sum(sub4);
sub5_tot = sum(sub5);
sub_tot_line = [sub1_tot,sub2_tot,sub3_tot,sub4_tot,sub5_tot];
depths = [dz(1),dz(1)+dz(2),dz(1)+dz(2)+dz(3),dz(1)+dz(2)+dz(3)+dz(4), dz(1)+dz(2)+dz(3)+dz(4)+dz(5)];

figure(1);

% make a plot of spectrally-resolved albedo:
plot(wvl,albedo,'linewidth',3);
xlabel('Wavelength (\mum)','fontsize',20);
ylabel('Albedo','fontsize',20);
set(gca,'xtick',0:0.1:5,'fontsize',16);
set(gca,'ytick',0:0.1:1.0,'fontsize',16);
xlim([0.3 2.5])
ylim([0,1])
grid on;
hold on;

% plot subsurface light field
figure(2);
plot(wvl, sub1, 'DisplayName','0 - 1 cm');
plot(wvl, sub2','DisplayName','1 - 3 cm');
plot(wvl, sub3,'DisplayName','3 - 5 cm');
plot(wvl, sub4,'DisplayName','5 - 7 cm');
plot(wvl, sub5,'DisplayName','7 - 9 cm');
xlabel('Wavelength (\mum)','fontsize',20);
ylabel('Planar Intensity Wm-2','fontsize',20);
xlim([0.3 2.5]);
ylim([0 0.03]);
set(gca,'xtick',0:0.1:5,'fontsize',16);
set(gca,'ytick',0:0.01:0.1,'fontsize',16);
legend;
hold on

%plot total energy intensity against depth
figure(3);
plot(depths,sub_tot_line);
xlabel('Depth beneath surface (m)','fontsize',20);
ylabel('planar intensity (Wm-2)','fontsize',20);
set(gca,'xtick',0:0.01:0.1,'fontsize',16);
set(gca,'ytick',0:0.1:1,'fontsize',16);
hold on

% report slope between blue and green reflectance
bgslope = ((albedo(26)-albedo(18))/0.075);

%Report albedo
alb_slr; % albedo over solar spectrum
%flx_abs_snw % absorbed energy in the snowpack
%heat_rt; % radiative heating rate in K/hr

snow_depth = sum(dz); % depth of snowpack incorporating all layers
temp_grad = ((heat_rt(1) - heat_rt(end)))/snow_depth; % temperature gradient through snowpack


end
