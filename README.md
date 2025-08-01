# PLAnET Microbial Model Description

The PLAnET model describes both phyllospheric population dynamics and net fluxes of microorganisms on the basis of simple meteorological variables.

SYNTAX:

`outstruc = planet(data,ustar_flag,depo_flag,params);`

PROGRAM STRUCTURE:

i) INPUTS:

The model needs at least three input arguments. A data matrix (data), a flag indicating if wind speed or friction velocity are included in the data matrix (ustar_flag) and a flag indicating if depositional effects are to be computed or not (depo_flag). A fourth argument (params)can be provided by the user if there's the need to use different model parameters from the ones derived by its original formulation. In the latter case, params is a vector of 11 values specifying: 1) Minimum Growth Temperature (°C) 
2) Maxmimum Growth Temperature (°C) 
3) Optimal Growth Temperature (°C) 
4) Lower population boundary (CFU m-2) 
5) Upper population boundary (CFU m-2) 6-8) three coefficients for the Gompertz equation governing upward fluxes 9) A scaling coefficient multiplying the growth of microorganisms 10-11) Slope and offset of the function relating LAI and the background airborne microbial concentration

The input data matrix has a specific format: it has a number of rows corresponding to the number of observations and four columns each one corresponding to a driving variable. The four columns are: 1) Temperature (°C) 
2) Pressure (Pa) 
3) Wind Speed / u* (m/s) 
4) Leaf area index 
5) Wind Speed (if only u* is provided in column 3, see below).

The ustar flag (ustar_flag) governs how the data matrix is treated by the model. If the flag == 1, it assumes that the third column in the matrix is already representing a friction velocity; if the flag == 0 it assumes that, instead, wind speed is given. The model converts wind speed in friction velocity on the basis of the height of the wind measurement in meter (z, in m) and the roughness length of the surface (z0, in m). These parameters are fixed at 10 and 0.15 m respectively but may be modified in the "MODEL CONSTANTS" section of the code by the user if needed. If ustar is directly given, the model needs wind speed as well to compute deposition. Wind speed should be added as a fifth column to the input matrix if deposition wants to be used.

The third input of the model is the flag governing the deposition. When depo_flag==1 the model calculates settling velocity and depositional fluxes using a diameter fixed at 3.3 microns. If the user wants to change it, it is possible to modify "dia" in the "MODEL CONSTANTS" section. When depo_flag==0 no depositional fluxes are computed and the net flux is never corrected. Besides gravitational settling velocity the model computes also the deposition due to interception/impaction (see function "slinnmod" at lines 325ff). General impaction/interception parameters can be tweaked under the function's "model parameters" (lines 336-340)

ii) OUTPUTS:

To avoid having the user input a dozen of output arguments, the model exits with only a single argument (outstruc). The latter is a structure containing all the main model variables: 1) outstruc.population = the microbial population in the phyllosphere (CFU m-2) 
2) outstruc.growth = the microbial growth rate 
3) outstruc.kmax = the variation in population cap due to leaf senescence 
4) outstruc.gross_out = gross outward flux (CFU m-2 s-1) 
5) outstruc.gross_in = gross inward flux (CFU m-2 s-1) 
6) outstruc.conc = background airborne concentration (CFU m-3) 
7) outstruc.net_flux = net microbial flux (CFU m-2 s-1) 
8) outstruc.vg = computed settling velocity due to gravity (m s-1) 
9) outstruc.vimp = computed settling velocity due to interception /impaction (m s-1) 
10) oustruc.vdep = overall settling velocity (vg + vimp) (m s-1) 
11) outstruc.dflag = deposition flag. If it's equal to 1 deposition is active at the given timestep. 
12) outstruc.dieout = dieout flow (CFU m-2)

iii) GENERAL GUIDELINES:

1) The model assumes an half-hourly timestep in the input, and that's why it divides the output values in fluxes by d=1800 in order to output fluxes in CFU m-2 s-1. If the user has other needs, the value of d in "MODEL CONSTANTS".
2) The model assumes to be starting from a moment when the microbial population is at its minimum, so it's advisable to have a dataset starting in winter months when LAI is at its minimum.
