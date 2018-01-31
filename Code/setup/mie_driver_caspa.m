% Driver for miecoated.m, originally written by Christian matzler (see
% Matzler, 2002).
% 
% inputs are:
%   rice : radius of inner ice sphere
%   rwater : radius of outer water sphere (from centre of ice sphere to
%           edge of water sphere)
% 
% This script was written to interface with the caspa setup scripts such
% that the variables rad and wat_rad are pulled in from
% snicar_melting_setup.m which creates snicar input variables for a melting
% snowpack driven by a grain size evolution code and an algal growth model.
% The resulting grain sizes and water coating thickness are parsed to this
% script which generates the relevant optical properties for the grains and
% water coatings generated and adds them to the working directory to be
% called from snicar during a CASPA run.

% For this reason, the run this script independently of CASPA, the
% variables rad and wat_rad are the inputs which are then aliased to rice
% and water for the call to mie_coated.m


% Outputs: 
% miecoated.m returns the following efficiency factors: [qext qsca qabs qb asy qratio]
% the driver convolves these with the particle dimensions to return the
% cross sections for extinction, scattering and absorption plus the
% asymmetry parameter, q ratio and single scattering albedo.
%
% Note that the call to miecoated.m includes a term 'opt', which can be set
% to the values 1,2 or 3. These are alternative methods for calculating the
% absorption efficiency which can perform differently in some rare cases.
% The default is set to 3 here as this seems to perform better for large
% size parameters, which are likely to be used when modelling melting ice.
%
% Also note that the code includes an interpolation regime. This is because
% the original code produced NaNs for a few wavelengths at certain size
% parameters, particularly in the mid NIR wavelengths.

% The code automatically saves new netcdfs files and adds them to the
% working directory


% Written by Joseph Cook, Feb 2018, University of Sheffield, UK
function [] = mie_driver_caspa(rad,wat_rad)

% set up variables
extinction=0;
scattering = 0;
absorption = 0;
asymmetry = 0;
q_ratio = 0;
ssa = 0;

% Define relevant parameters: Only rice and rwater are user defined.
rice = rad; % inner sphere diameter in um
rwater = rad+wat_rad; % outer sphere diameter in um (i.e. total coated sphere, not water layer thickness)

% calculate volume and density of sphere

XSArea_inner = pi*((rice)^2); % cross sectional area of ice core
XSArea_outer = pi*((rwater)^2) - XSArea_inner; % cross sectional area of water layer sphere
TotalXS = pi*((rwater)^2);


WatDensity = 999; % density of water at 1 degree C in kg m-3
IceDensity = 934; % density of ice in kg m-3

IceVol = 4/3 * pi * (rice)^3;
WatVol = (4/3 * pi * (rwater)^3) - IceVol;
TotalVol = 4/3 * pi * (rwater)^3;


IceMass = IceVol * IceDensity;
WatMass = WatVol * WatDensity;
TotalMass = IceMass + WatMass;


%define wavelength range
WL= 0.305:0.01:5;


for i = (1:1:470)
    % size parameters for inner and outer spheres
   x = 2 * pi * rice / (WL(i));

   y = 2 * pi * rwater / (WL(i));
   
    % read in refractive indices from refracICE and refracWATER per wavelength
   [RN CN] = refracICE(WL(i));
   m1 = complex(RN,CN);
   [RN2 CN2] = refracWATER(WL(i));
   m2 = complex(RN2,CN2);
   % call Miecoated and return efficiencies
   [qext qsca qabs qb asy qratio] = Miecoated(m1,m2,x,y,3);
   %append efficiencies to lists
   extinction(i)=qext;
   scattering(i) = qsca;
   absorption(i) = qabs;
   backscattering(i) = qb;
   asymmetry(i) = asy;
   q_ratio(i) = qratio;
   ssa(i) = qsca/qext;
end

% replace any possible nans with values estimated by cubic interpolation
extinction = fillmissing(extinction,'spline')
scattering= fillmissing(scattering,'spline')
absorption= fillmissing(absorption,'spline')
backscattering= fillmissing(backscattering,'spline')
asymmetry= fillmissing(asymmetry,'spline')
q_ratio= fillmissing(q_ratio,'spline')
ssa= fillmissing(ssa,'spline')

% Instability exists for long wavelengths and large grain sizes.
% Experimentation showed this could occur for coated grains > 1300 microns.
% The following loop protects against unrealistic values at those long
% wavelengths by interpolating from the rest of the data for the final 70
% wavelengths (i.e. from 4.3 - 5 um)

if rwater > 1300
    
    iWL = WL(1:400);
    iextinction = extinction(1:400);
    i_extinction = interp1(iWL,iextinction,[4.3:0.01:5],'linear','extrap');
    extinction(400:470) = iextinction(400:end)
    plot(extinction)

    iscattering = scattering(1:400);
    i_scattering = interp1(iWL,iscattering,[4.3:0.01:5],'linear','extrap');
    scattering(400:470) = iscattering(400:end)
    figure
    plot(scattering)
    
    iabsorption = absorption(1:400);
    i_absorption = interp1(iWL,iabsorption,[4.3:0.01:5],'linear','extrap');
    absorption(400:470) = iabsorption(400:end)
    figure
    plot(absorption)
    
    iSSA = ssa(1:400);
    i_SSA = interp1(iWL,iSSA,[4.3:0.01:5],'linear','extrap');
    ssa(400:470) = iSSA(400:end)
    figure
    plot(ssa)
    
    iASM = asymmetry(1:400);
    i_ASM = interp1(iWL,iASM,[4.3:0.01:5],'linear','extrap');
    asymmetry(400:470) = iASM(400:end)
    figure
    plot(asymmetry)
    
end


% calculate cross sections from efficiency factors

ExtXC = (extinction .* TotalXS);
ScaXC = (scattering .* TotalXS);
AbsXC = (absorption .* TotalXS);

ExtXCvol = (extinction.*TotalVol);
ScaXCvol = (scattering.*TotalVol);
AbsXCvol = (absorption.*TotalVol);

ExtXCmass = (extinction.*TotalMass);
ScaXCmass = (scattering.*TotalMass);
AbsXCmass = (absorption.*TotalMass);

% print fraction of sphere made up of water by mass and volume
water_frac_mss = 100-(IceMass/TotalMass)*100;
water_frac_vol = 100-(IceVol/TotalVol)*100;
ice_frac_mss = 100 - water_frac_mss;
ice_frac_vol = 100 - water_frac_vol;


part_dens = ((IceDensity*ice_frac_mss/100)+(WatDensity*water_frac_mss/100)); % density of particle as average weighted by mass of components


%%%%%%%%%%%%%%%%%%%%% CREATE AND POPULATE NETCDF  %%%%%%%%%%%%%%%%%%

% Copy existing netcdf to ensure formatting is consistent

filename = sprintf('ice_wrn_coated_%d_%d.nc',rad,wat_rad);

copyfile ('ice_wrn_coated_1000_100.nc', filename)

rr(1:200) = rad*10e-6 + wat_rad*10e-6; %radius in m

ncwrite(filename,'rds',rr) % particle radius
ncwrite(filename,'ext_xsc',ExtXC) % extinction cross section
ncwrite(filename,'sca_xsc',ScaXC) % scattering cross section
ncwrite(filename,'abs_xsc',AbsXC) % absorption cross section
ncwrite(filename,'ext_cff_mss',ExtXCmass) % mass extinction cross section
ncwrite(filename,'sca_cff_mss',ScaXCmass) % mass scattering cross section
ncwrite(filename,'abs_cff_mss',AbsXCmass) % mass absorption cross section
ncwrite(filename,'ext_cff_vlm',ExtXCvol) % volume extinction cross section
ncwrite(filename,'sca_cff_vlm',ScaXCvol) % volume scattering cross section
ncwrite(filename,'abs_cff_vlm',AbsXCvol) % volume absorption cross section
ncwrite(filename,'ss_alb',ssa) % single scattering albedo
ncwrite(filename,'asm_prm',asymmetry) % assymetry parameter

% PSD variables are assigned but are not used by SNICAR (only uses Reff)
% This is a legacy thing from the original SNICAR code that included some
% PSD diagnostics. Particle density is calculated and assigned. 

ncwrite(filename,'rds_swa',rad*10e-6 + wat_rad*10e-6) % surface weighted radius (analytic)
ncwrite(filename,'rds_swr',rad*10e-6 + wat_rad*10e-6) % surface weighted radius (resolved)
ncwrite(filename,'rds_nma',1e-4) % analytic number-median radius
ncwrite(filename,'gsd',1.5) % geometric SD of lognormal distribution
ncwrite(filename,'prt_dns',part_dens) % particle density

end


