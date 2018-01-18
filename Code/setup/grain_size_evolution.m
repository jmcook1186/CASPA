 
% We've altered the functional form of the best-fit equations from our 2006 JGR paper, so they now represent the evolution of effective radius 
% instead of specific surface area (SSA).  The primary reason for doing so was that small errors in low SSA values produce large differences in 
% effective radius (since they are inversely related), and hence optical properties.  The parametric equation is also now formulated in terms of 
% the rate of change of effective radius, to accommodate evolving snow states:
% 
% dr/dtime = (drdt0*(tau/((r-r0)+tau))^(1/kappa))
% 
% where drdt0, tau, and kappa are retrieved from the look-up tables, depending on snow temperature, temperature gradient, and density. "r" is the
% current effective radius in microns, "r0" is the fresh snow effective radius (so that r-r0 is positive definite), and drdt0 is the initial rate 
% of change of effective radius (in microns per hour).
% 
% In these files, drdt0, tau, and kappa all have dimensionality: [T dTdz sno_dns], where T is temperature in Kelvins, dTdz is temperature gradient (K/m), 
% and sno_dns is snow density in kg/m3.  The ranges and values of these three dimensions (corresponding to the 11x31x8 parameter matrices) are defined in 
% 1-D arrays (T, dTdz, and sno_dns).
% 
% There are three files in each format, containing parameters for the evolution of snow with initial SSA values of 60, 80, and 100 m2/kg. 
% In the NCAR Community Land Model we apply SSA0=60 parameters, but other literature suggests that fresh snow generally has higher SSA than this. 
% Conversion between effective radius and SSA is simply: r=3/(SSA*917).
% 
% Finally, in CLM we also include crude representations of aging associated with liquid water content and re-freezing of liquid water.  
% Basically these partial aging terms are added to the dry aging rate that is calculated from the lookup table.  
% The full scheme is described in the CLM technical document, starting on p.38: http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf


% load file and set variables (later will be pulled from snicar)
function [r_new] = grain_size_evolution(temp,TG,density,r0,r,doubling_time,fliq,f_ref)

load('drdt_bst_fit_60.mat');

temp = temp;% temperature in K
TG = TG; % temp gradient
density = density; % snow density in kgm-3
r0 = r0;
r = r; % ice crystal diameter in microns
doubling_time = doubling_time;
fliq = fliq;
f_ref = f_ref; % fraction of refrozen ice, assumed reff = 1500 um
C = 4.22e-13; % constant (Flanner and Zender 2006)

% find index of closest match for each variable in database

[error_sno_dns, idx_sno_dns] = min(abs(sno_dns - density(1)));
[error_T, idx_T] = min(abs(T - temp(1)));
[error_dTdz, idx_dTdz] = min(abs(dTdz - TG(1)));

% index into each LUT to retrieve value for kappa, tau and drdt0

kappa_v = kappa(idx_T,idx_dTdz,idx_sno_dns);
tau_v = tau(idx_T,idx_dTdz,idx_sno_dns);
drdt0_v = drdt0(idx_T,idx_dTdz,idx_sno_dns);

dry_aging = (drdt0_v*(tau_v/((r-r0)+tau_v))^(1/kappa_v));
wet_aging = (10e18 * C * fliq^3) / (4* pi * r^2);
refrozen = 1500 * f_ref;

dr_dtime = dry_aging + wet_aging + f_ref; % rate of change in microns/hr

dr_day = dr_dtime * 24; %change over one day
dr_timestep = dr_day * doubling_time; % change over timestep

r_new = r + dr_timestep;




