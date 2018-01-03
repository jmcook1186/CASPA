# CASPA : Cellular Automaton for SnowPack Albedo

README written and maintained by J Cook, University of Sheffield

Initial commits contain 4 versions of CASPA code previously held locally by J Cook. These four versions were produced before version control migrated to Github, so updates are marked in comments in scripts.

The code is written in Matlab and is designed to provide a framework for modelling snow albedo using radiative transfer by invoking the BioSNICAR model (Cook et al., 2017) which was an adapted form of Mark Flanner's original SNICAR model enabling the incorporation of biological impurities. Here, a cellular automaton is developed that uses BioSNICAR to calculate the albedo of a snowpack in each cell in a grid i x j, and evolves the snowpack according to biological growth rules dervied from empirical experiments and existing literature. 

This code is in active development and is not yet published. There are no permissions granted for copying, applying or using this code.


# OVERVIEW

This model simulates a snow surface with physical/optical properties
defined by the user. This surface is represented as a cellular automaton,
i.e. an x,y grid composed of discrete cells that update based upon certain
functions as the model advances through time.
There are various modes of operation. The simplest is to simulate an
algal bloom that grows from an initial configuration defined by the user.
The user defines the % coverage at t=0. At each timestep, algal growth is
simulated using a probabalistic function, where growth in situ has a
given likelihood, spread of the algae into neighbouring cells has a given
likelihood, and spontaneous initiation of a new bloom has a give
likelihood. Growth is represented by increasing the cell value - i.e. 0
represents clean ice, 1-10 represent increasing mixing ratio of algae in
the upper surface. A carrying capacity is simulated wherein the cell
value cannot exceed a value of 10. Once the cell value reaches 10 the
bloom can only grow by spreading. If all the neighbouring cells are also
at carrying capacity, the algae will not grow.

At each timestep, the grid is thereby updated. The value of each cell is
associated with a call to a specific instance of the radiative transfer
scheme BioSNICAR. These instances are pre-coded as driver routines
saved in the working directory. In this default version, each instance is
identical except for the mixing ratio of algae in the upper 3mm of ice.
The broadband and spectral albedo of each pixel in the grid is thereby 
modelled using BioSNICAR and a spatially-integrated albedo computed 
using a grid-mean. The albedo is used to calculate the radiative forcing
following Ganey et al (2017: Nat Geoscience) and Dial et al (2018: FEMS)
by multipling the incoming irradiance by 1-albedo for clean and 'mixed'
surfaces, the difference between which provides the RF caused by
impurities.

One timestep is dimensionless, except that successive SNICAR instances
have 2x algal biomass in upper layer, meaning 1 timestep is equal to the
doubling time of the algae, making the time dimension tractable. For a
doubling time of 3 days, one timestep = 3 days. Albedo change / doubling
time (days) = albedo change per day.
Also, in the default version the algal growth is uninterupted, so left
long enough the algae will have even coverage at carrying capacity.
However, there is the option to interrupt the growth, simulating a
rainfall event for example. Could potentially simulate nutrient
limitation by increasin or decreasing the doubling times for the blooms.
A nutrient nutrient pulse could thereby be simulated.

Possibility to model indirect effects: could increase grain size along
with algal biomass, according to field observations (measured biomass vs
depth of 1030nm absorption feature). 

# Model Development Prior to Initial Commit to Github

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

# TODOS
       
1. Deal with edge effects (-1s and +1s can push the indexing out of grid
range, currently counter starts at 2 and ends at i,j -1, but this
leads to a border that doesn't update)
    
2. plot(meangrid) after runs and analyse temporal evolution
     
3. Could add an additional randi() dice roll to randomly reduce value
of cells, simulating rainfall washaway. Alternatively, a sudden dust
deposition event...? These don't have to be random, they could start
on a particular ticker value (e.g. at 20th timestep: rain).
        
4. Can the initial grid be set using UAV por satellite imagery?
    
5. How should I deal with dust? constant background value? 
