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
function [] = mie_driver_caspa(rad)

% set up variables
extinction=0;
scattering = 0;
absorption = 0;
asymmetry = 0;
q_ratio = 0;
ssa = 0;

% Define relevant parameters: Only rice and rwater are user defined.
rice = rad; % inner sphere diameter in um

% calculate volume and density of sphere

XSArea = pi*((rice)^2); % cross sectional area of ice core

IceDensity = 934; % density of ice in kg m-3

IceVol = 4/3 * pi * (rice)^3;


IceMass = IceVol * IceDensity;


%define wavelength range
WL= 0.305:0.01:5;

for i = (1:1:470)
   [qext, qsca, qabs, qb, asy, qratio]=MieIce(WL(i),rice);
   extinction(i)=qext;
   scattering(i) = qsca;
   absorption(i) = qabs;
   backscattering(i) = qb;
   asymmetry(i) = asy;
   q_ratio(i) = qratio;
   ssa(i) = qsca/qext;
end

% calculate cross sections from efficiency factors

ExtXC = (extinction .* XSArea);
ScaXC = (scattering .* XSArea);
AbsXC = (absorption .* XSArea);

ExtXCvol = (extinction.*IceVol);
ScaXCvol = (scattering.*IceVol);
AbsXCvol = (absorption.*IceVol);

ExtXCmass = (extinction.*IceMass);
ScaXCmass = (scattering.*IceMass);
AbsXCmass = (absorption.*IceMass);


part_dens = 914; % density of ice

%%%%%%%%%%%%%%%%%%%%% CREATE AND POPULATE NETCDF  %%%%%%%%%%%%%%%%%%

% Copy existing netcdf to ensure formatting is consistent


rad4dig = sprintf('%04d',rad);
filename = sprintf('ice_wrn_%s.nc',rad4dig);
copyfile ('ice_wrn_1000.nc', filename)



rr(1:200) = rad*10e-6; %radius in m

ncwrite(filename,'rds',rr) % particle radius
ncwrite(filename,'ext_xsc',ExtXC) % extinction cross section
ncwrite(filename,'sca_xsc',ScaXC) % scattering cross section
ncwrite(filename,'abs_xsc',AbsXC) % absorption cross section
ncwrite(filename,'ext_cff_mss',ExtXCmass) % mass extinction cross section
ncwrite(filename,'sca_cff_mss',ScaXCmass) % mass scattering cross section
ncwrite(filename,'abs_cff_mss',AbsXCmass) % mass absorption cross section
ncwrite(filename,'ss_alb',ssa) % single scattering albedo
ncwrite(filename,'asm_prm',asymmetry) % assymetry parameter

% PSD variables are assigned but are not used by SNICAR (only uses Reff)
% This is a legacy thing from the original SNICAR code that included some
% PSD diagnostics. Particle density is calculated and assigned. 

ncwrite(filename,'rds_swa',rad*10e-6) % surface weighted radius (analytic)
ncwrite(filename,'rds_swr',rad*10e-6) % surface weighted radius (resolved)
ncwrite(filename,'rds_nma',1e-4) % analytic number-median radius
ncwrite(filename,'gsd',1.5) % geometric SD of lognormal distribution
ncwrite(filename,'prt_dns',part_dens) % particle density

end


