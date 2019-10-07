# TheDiaTo v1.0
A collection of diagnostic for the study of the thermodynamics of the climate system

## Author
Valerio Lembo
(Meteorological Institute, Hamburg University - valerio.lembo@uni-hamburg.de)

## Introduction

The tool consists of four modules; one for the
computation of energy budgets and transports, one for the water mass budgets
(and related meridional transports), one for the Lorenz Energy Cycle (LEC),
one for the material entropy production.

The first module is run by default, the others are optional. If the lsm option
is set to true, the module 1 and the module 2 will be run with additional
separate results over land and oceans.

- MODULE 1 (default):
Earth's energy budgets from radiative and heat fluxes at Top-of-Atmosphere,
at the surface and in the atmosphere (as a residual).
Meridional transports, magnitude and location of the peaks in each
hemisphere (only for heat transports) are also computed.
The baroclinic efficiency is computed from TOA energy budgets, emission
temperature (in turn retrieved from OLR) and near-surface temperature.

- MODULE 2 (optional):
Water mass and latent energy budgets and meridional transports are computed
from latent heat fluxes, snowfall and rainfall precipitation fluxes. Magnitude
and location of the peaks in each hemisphere (only for heat transports) are
also computed, as for module 1.

- MODULE 3 (optional):
The Lorenz Energy Cycle (LEC) is computed in spectral components from near-
surface temperatures, temperatures and the three components of velocities
over pressure levels.
The storage and conversion terms are directly computed, the sources and
sinks are retrieved as residuals.
Components are grouped into a zonal mean, stationary and transient eddy
part.

- MODULE 4 (optional):
The material entropy production is computed using the indirect method, the
direct method or both (following Lucarini et al., 2014).
For the indirect method a vertical and a horizontal component are provided.
For the direct method, all components are combined, related to the
hydrological cycle (attributable to evaporation, rainfall and snowfall
precipitation, phase changes and potential energy of the droplet), to the
sensible heat fluxes and to kinetic energy dissipation. For the latter the
LEC computation is required, given that the strength of the LEC can be
considered as equal to the kinetic energy dissipated to heating. If the option
for module 3 is set to false, a reference value for the material entropy
production related to the kinetic energy dissipation is provided.

## Prerequisites

- Python v3.6x;
- CDO bindings for Python v1.5.x (update to latest available version);
- Matplotlib v2.x (at least 2.2);
- Numpy v1.15.x;
- Netcdf4 v1.4.2; 

The program has been tested in an environment built through Anaconda3
repository: this might be necessary if the user is working on a server with
no admin permissions. The above mentioned packages are retrieved from
conda-forge collections. In case any question or problem arises, please contact
the developers.

## Content of the Package
- thermodyn_diagnostics.py: the main script;
- computations.py: contains all core computations;
- fluxogram.py: contains the instructions to produce a *.png picture with the
                diagram of the Lorenz Energy Cycle (LEC);
- fourier_coefficients.py: contains the script for producing Fourier
  coefficients from gridded data;
- lorenz_cycle.py: contains the functions for computation of the LEC;
- mkthe.py: contains the functions for computation of auxiliary variables;
- namelist.py: a namelist script to be adapted by the user;
- plot_script.py: contains functions for all plots (except the LEC diagram);

## Usage

1. Obtain the datasets and store them as follows:

    a. The following variables are needed at the monthly resolution or higher
       (required variable name in brackets):
     - TOA shortwave radiation downwards (rsdt);
     - TOA shortwave radiation upwards (rsut);
     - TOA longwave radiation upwards (rlut);
     - Surface shortwave radiation downwards (rsds);
     - Surface shortwave radiation upwards (rsus);
     - Surface longwave radiation downwards (rlds);
     - Surface longwave radiation upwards (rlus);
     - Surface upward turbulent latent heat fluxes (hfls);
     - Surface upward turbulent sensible heat fluxes (hfss);
     - Precipitation flux (pr);
     - Snowfall flux (prsn);
     - Surface air pressure (ps);
     - Surface temperature (ts);
     - Specific humidity on pressure levels (hus);
     
     N.B. Precipitation fluxes are expected to have kg m-2 s-1 as u.o.m.
     
    b. The following variables are needed at the daily resolution or higher:
     - Near-surface temperature (tas);
     - Near-surface (or 10m) zonal velocity (uas);
     - Near-surface (or 10m) meridional velocity (vas);
     - Air temperature (on pressure levels) (ta);
     - Horizontal velocity (on pressure levels) (ua);
     - Meridional velocity (on pressure levels) (va);
     - Vertical velocity (on pressure levels) (wap);
     
    c. The following dimensions are accepted (dimension name in brackets):
     - Longitude (lon);
     - Latitude (lat);
     - Pressure levels (plev);
     - Time (time);
     
       N.B. The latitude must be S->N (the CDO command 'invertlat' may help adapting
       data to this requirement). The datasets must be global, i.e. spanning the
       whole Earth.
       
       N.B. A calendar time must be defined and the time axis must refer to that. You may
       use the combination of 'setreftime' and 'settaxis' CDO commands in order to prepare
       the datasets. In doing so, the time resolution shall be correctly specifiend (see the
       CDO manual for more information).
       
       N.B. The grids are expected to be recognised structured lonxlat grids (either regular or gaussian).
       If this is not the case, you may want to remap to recognised grids, using either remapbil or 
       remapcon CDO commands (see the CDO manual for specific instructions).
     
    d. In case variable and dimension names do not comply to these requirements
       you may use the CDO command (for variables): 
       'chname,namein,nameout filein.nc fileout.nc' 
       and the NCO command (for dimensions): 
       'ncrename -h -d namein,nameout -v .namein,nameout filename.nc'
    
    e. Create a folder in the input directory for each model you need to
       analyse, named after the model.
    
    f. Provide each field in a different file. The file name must start with the name of the variable.
    
    g. If the option 'lsm' is set to 'True', a land-sea mask shall be provided as a land area fraction (from 0 to 100) in a          separate folder, named as 'ls_mask', at the same level as the folder containing the model outputs for a single model.
       The land-sea mask is expected to be a NetCDF file whose name format is 'lsmask_<model_name>.nc'.
       Example:
        
        model_data/
        model_data/model_name/*.nc
        model_data/ls_mask/lsmask_<model_name>.nc

2: The companion script 'namelist.py' is available for the user in order to
   specify the directory paths. The following options can also be set:
   - wat: if set to true, the program will compute the water mass and
          latent energy budget,
   - lec: if set to true, the program will compute the Lorenz Energy Cycle
              (LEC) averaged on each year;
   - entr: if set to true, the program will compute the material entropy
           production (MEP);
   - met: if set to 1, the program will compute the MEP with the indirect
          method, if set to 2 with the direct method, if set to 3, both
          methods will be computed and compared with each other;
              
4: Run the tool by typing: 'python thermodyn_diagnostics.py'. Global metrics
   and additional information is printed out in the log file indicated by the
   user in the 'namelist.py' script.

## Output

The output directory contains the following NetCDF files:

- (output directory):
   - atmos_transp_mean_<model_name>.nc
   - latent_transp_mean_<model_name>.nc
   - ocean_transp_mean_<model_name>.nc
   - total_transp_mean_<model_name>.nc
   - wmb_transp_mean_<model_name>.nc

   contain annual mean meridional sections of heat transports in the
   atmosphere, oceans, and as a total; latent energy transports and water
   mass transports;

- (output directory)/<model_name>:
  - <model_name>_atmb.nc
  - (<model_name>_latent.nc; if wat is set to true)
  - <model_name>_surb.nc
  - <model_name>_toab.nc
  - (<model_name>_wmb.nc; is wat is set to true)
  
  contain annual mean 2D fields of energy budget, latent heat and water
  mass budgets;

  - <model_name>_barocEff.nc
  
  contains the evolution of annual mean baroclinic efficiency (Lucarini et al., 2011).

  (if entr is set to true):
  - <model_name>_evap_entr.nc (if met is set to 2 or 3)
  - <model_name>_horizEntropy.nc (if met is set to 1 or 3)
  - <model_name>_pot_drop_entr.nc (if met is set to 2 or 3)
  - <model_name>_rain_entr.nc (if met is set to 2 or 3)
  - <model_name>_sens_entr.nc (if met is set to 2 or 3)
  - <model_name>_snow_entr.nc (if met is set to 2 or 3)
  - <model_name>_snowmelt_entr.nc (if met is set to 2 or 3)
  - <model_name>_verticalEntropy.nc (if met is set to 1 or 3)
  
  contain the evolution of annual mean components of the material entropy production.
  
- (plots directory):
  - meridional_transp.png: contains the model inter-comparison of total, atmospheric and oceanic of meridional
                           sections in zonally average meridional heat transports;
  - scatters_summary.png: contains the scatter plots of model intercomparisons of various metrics retrieved in
                          the program;
  - scattes_variability.png: contains scatter plots of model intercomparisons between TOA, atmospheric and surface
                             global mean energy budgets and their inter-annual variability;

- (plots directory)/<model_name>:
  - <model_name>_atmb_timeser.png: the atmospheric budget annual mean global and hemispheric time series;
  - <model_name>_energy_climap.png: the TOA, atmospheric and surface climatological mean fields;
  - <model_name>_heat_transp.png: the meridional sections of total, atmospheric and oceanic meridional
                                  heat transports (implied from energy budgets);
  - <model_name>_latent_climap.png: the climatological mean latent heat field;
  - <model_name>_latent_timeser.png: the latent heat annual mean global and hemispheric evolutions;
  - <model_name>_latent_transp.png: the meridional section of annual mean meridional latent heat transport;
  - <model_name>_scatpeak.png: the scatter plots of atmospheric vs. oceanic peak magnitude in both hemispheres;
  - <model_name>_sevap_climap.png: the annual mean field of material entropy production due to evaporation;
  - <model_name>_smelt_climap.png: the annual mean field of material entropy production due to snow melting;
  - <model-name>_spotp_climap.png: the annual mean field of material entropy production due to potential energy
                                   of the droplet;
  - <model_name>_srain_climap.png: the annual mean field of material entropy production due to rainfall precipitation;
  - <model_name>_ssens_climap.png: the annual mean field of material entropy production due to sensible heat fluxes;
  - <model_name>_ssnow_climap.png: the annual mean field of material entropy production due to snowfall precipitation;
  - <model_name>_surb_timeser.png: the surface budget annual mean global and hemispheric time series;
  - <model_name>_sver_climap.png: the annual mean field of vertical material entropy production through the indirect method;
  - <model_name>_toab_timeser.png: the TOA budget annual mean global and hemispheric time series;
  - <model_name>_wmb_climap.png: the climatological mean water mass budget field;
  - <model_name>_wmb_timeser.png: the water mass annual mean global and hemispheric evolutions;
  - <model_name>_wmb_transp.png: the meridional section of annual mean meridional water mass transport;

- (plots directory)/<model_name>/LEC_results:
  - <model_name>_<year>_lec_diagram.png: the flux diagram for the annual mean LEC cycle in a specific year;
  - <model_name>_<year>_lec_table.txt: the table containing the storage and conversion terms for the annual mean
                                       LEC cycle in a specific year;

## References

Lembo, V., Lunkeit, F., and Lucarini, V.: TheDiaTo (v1.0) – A new diagnostic tool for water, energy and entropy budgets in climate models, Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2019-37, in review, 2019

## Versions
2019-04-29: a stand-alone version of TheDiaTo v1.0 is branched from ESMValTool v2.0b repository;
