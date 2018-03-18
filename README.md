# CASPA : Cellular Automaton for SnowPack Albedo

README written and maintained by J Cook, University of Sheffield

The code is written in Matlab and is designed to model snow albedo using radiative transfer by invoking the BioSNICAR model (Cook et al., 2017) in a cellular automaton framework. BioSNICAR is used to calculate the albedo of a snowpack in each cell in a grid i x j. At each timestep the algal biomass is updated according to a user defined doubling time and a grain size evolution model adapoted from the Community Land Model and Flanner and Zender(2006) is used to update the snow grain sizes and layer thicknesses, driven by energy fluxes predicted by te radiative transfer model.

This code is in active development and is not yet published. There are no permissions granted for copying, applying or using this code.


# OVERVIEW

This model simulates a snow surface with initial physical/optical properties
defined by the user. This surface is represented as a cellular automaton,
i.e. an x,y grid composed of discrete cells that update based upon certain
functions as the model advances through time.

The model simulates an algal bloom that grows from an initial configuration
defined by the user. The user defines the % coverage at t=0. At each timestep,
algal growth is simulated using a probabalistic function, where growth in situ
has a given likelihood and spread of the algae into neighbouring cells
(Moore neighbourhood) has a given likelihood. Growth is represented by 
increasing the cell value - i.e. 0 represents clean ice, 1-10 represent 
increasing mixing ratio of algae in the upper surface. A carrying capacity 
is simulated wherein the cell value cannot exceed a value of 10. Once the 
cell value reaches 10 the bloom can only grow by spreading. If all the 
neighbouring cells are also at carrying capacity, the algae will not grow.

At each timestep, the grid is thereby updated. The value of each cell is
associated with a call to a specific instance of the radiative transfer
scheme BioSNICAR. BioSNICAR is run offline coupled to a grain size evolution
model that uses the energy fluxes in each vertical layer of the snowpack
predicted by BiOSNICAR. The grain size evolution model accounts for wet and 
dry grain growth and accumulation of liquid water. Liquid water is 
generated whenever there is excess energy that raises the snow temperature
abive 273.15K. This water can pecolate into lower layers. If the lower layers
have sub-freezing temperatures the water will refreeze. refrozen water is
assumed to have a radius of 1500 microns. Liquid water that remains in situ
is modelled by increasing the optical radius of the ice grains. The optical
properties of the snow is therbey updated for each instance of BioSNICAR and
the broadband and spectal albedo are returned for each surface class. These
are added to a lookup table accessible by CASPA as it updates the cellular
model each timestep. Running CASPA with the grain evolution model ON and OFF 
and differencing the twp outputs provides a measure of direct and indirect 
effects of algae on albedo and radiative forcing.

A spatially-integrated spectral and broadband albedo for the entire grid is 
computed using a 2D-mean. Biomass is calculated from the volume, density and 
mixing ratio of impurities in the upper snow layer.

The spectral albedo and spetcral irradiance is used to calculate the radiative
forcing following Ganey et al (2017: Nat Geoscience) and Dial et al (2018: FEMS)
by multipling the incoming irradiance by 1-albedo for clean and 'mixed'
surfaces, the difference between which provides the RF caused by
impurities.

One timestep is dimensionless, except that successive SNICAR instances
have 2x algal biomass in upper layer, meaning 1 timestep is equal to the
doubling time of the algae, making the time dimension tractable. For a
doubling time of 3 days, one timestep = 3 days. Albedo change / doubling
time (days) = albedo change per day.


# RUNNING THE MODEL

Running the model in default mode is highly automated and controlled from a
driver script (CASPA_driver.m).

The CASPA repository should be cloned or downloaded and all files
saved to the working directory. 

The user-defined variables are assigned in the driver. The grain size evolution
code can be turned on or off and the properties of the snowpack defined.

The scavenging rate for mineral dusts is user defined. To have a dust free
model run, set the initial mixing ratio (y) to zero. To have a constant dust
mixing ratio throughout the run, set the scavenging rate to 1.

Running the driver first calls the setup functions. This will populate the 
workspace with all relevant datasets. 
Then CASPA is called (CASPA_V01.m at time of writing README).

The script will produce plots of the snow surface showing the growing algal 
blooms, albedo and IRF against time, biomass and coverae against time, and
biomass and albedo against time.

To tweak the inputs, user-defined variables can be found all together near 
the start of the script. For more in-depth changes, such as SNICAR conditions,
bio-optical parameters etc, refer to the main documentation provided in this
repository.


# Model Development Prior to Initial Commit to Github
These notes detail the state of the code prior to the first upload to GitHub,
for reference.

## CASPA v 0.01

This version is the cellular automaton only, i.e it evolves a surface using
classifiers but the model does not use these classifiers to call SNICAR. This
was the very first stage of model development.


## CASPA v 0.02

This version was the first to call SNICAR from the cellular model. It looped
through each cell in the grid and called SNICAR each time. This was slow and was
raised as an urgent TODO for v0.03. 


## CASPA v 0.03

This version no longer runs SNICAR per cell. Instead, all the instances of SNICAR used
in the model are run once at the beginning and the spectral and broadband
data appended to lists. These are then used as a lookup table for the
cellular model. Of course, this has drastically reduced the run time.

Bug fix in final loop - previously able to increase cell values above
10. This was because when an individual cell reached carrying capacity
the algae grew by spreading but was able to spread into cells whose
values > 10. This is now fixed so no cell can exceed the carrying
capacity.

Condense loops - now single loop for if rand>50 OR cell value 
 =10 (previously two sequential loops)


## CASPA v 0.04

Updated colormap to show clean snow as white, then grades of red as snow algae
bloom. This version is the first to include radiative forcing calculations in
addition to albedo. This is achieved by differencing the spectral albedo of the
simulated snowpack and a hypothetical clean snowpack then multiplying
1 - albedo by incoming spectral irradiance and normalising in space/time. There were
also some minor changes to variable names etc. to imporve the readability.
