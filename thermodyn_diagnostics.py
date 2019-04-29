r"""MAIN PROGRAM.

TheDiaTo - The diagnostic tool for climate system thermodynamics.

Author
Valerio Lembo
(Meteorological Institute, Hamburg University - valerio.lembo@uni-hamburg.de)

Contributors
Frank   Lunkeit
(Meteorological Insitute, Hamburg University - f.lunkeit@uni-hamburg.de)
Nikolay Koldunov
(MARUM/AWI, nikolay.koldunov@awi.de, Germany)

Project
CRC - TRR 181 "Energy transfers in Atmosphere and Ocean"

#############################################################################

SOFTWARE DESCRIPTION

The tool consists of three modules; one for the
computation of energy budgets and transports, one for the water mass budgets
(and related meridional transports), one for the Lorenz Energy Cycle (LEC),
one for the material entropy production.

The first module is run by default, the others are optional. If the lsm option
is set to true, the module 1 and the module 2 will be run with additional
separate results over land and oceans.

- MODULE 1 (default)
Earth's energy budgets from radiative and heat fluxes at Top-of-Atmosphere,
at the surface and in the atmosphere (as a residual).
Meridional transports, magnitude and location of the peaks in each
hemisphere (only for heat transports) are also computed.
The baroclinic efficiency is computed from TOA energy budgets, emission
temperature (in turn retrieved from OLR) and near-surface temperature.

- MODULE 2 (optional)
Water mass and latent energy budgets and meridional transports are computed
from latent heat fluxes, snowfall and rainfall precipitation fluxes. Magnitude
and location of the peaks in each hemisphere (only for heat transports) are
also computed, as for module 1.

- MODULE 3 (optional)
The Lorenz Energy Cycle (LEC) is computed in spectral components from near-
surface temperatures, temperatures and the three components of velocities
over pressure levels.
The storage and conversion terms are directly computed, the sources and
sinks are retrieved as residuals.
Components are grouped into a zonal mean, stationary and transient eddy
part.

- MODULE 4 (optional)
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

PREREQUISITES

Python v3.6x
CDO bindings for Python v1.5.x (update to latest available version)
Matplotlib v2.x (at least 2.2)
Numpy v1.15.x
Netcdf4 v1.4.2 

The program has been tested in an environment built through Anaconda3
repository: this might be necessary if the user is working on a server with
no admin permissions. The above mentioned packages are retrieved from
conda-forge collections. In case any question or problem arises, please contact
the developers.

CONTENT OF THE PACKAGE
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

USAGE

1: Obtain the datasets and store them as follows:
    a. The following variables are needed at the monthly resolution or higher
       (required variable name in brackets):
     - TOA shortwave radiation downwards (rsdt);
     - TOA shortwave radiation upwards (rsut);
     - TOA longwave radiation upwards (rlut);
     - Surface shortwave radiation downwards (rsds);
     - Surface shortwave radiation upwards (rsus);
     - Surface longwave radiation downwards (rlds);
     - Surface longwave radiation upwards (rlus);
     - Surface turbulent latent heat fluxes (hfls);
     - Surface turbulent sensible heat fluxes (hfss);
     - Surface temperature (tas);
     - Specific humidity on pressure levels (hus);
    b. The following variables are needed at the daily resolution or higher:
     - Near-surface temperature;
     - Near-surface (or 10m) zonal velocity;
     - Near-surface (or 10m) meridional velocity;
     - Air temperature (on pressure levels);
     - Horizontal velocity (on pressure levels);
     - Meridional velocity (on pressure levels);
     - Vertical velocity (on pressure levels);
    c. The following dimensions are accepted (dimension name in brackets):
     - Longitude (lon);
     - Latitude (lat);
     - Pressure levels (plev);
     - Time (time);
     The latitude must be N->S (the CDO command 'invertlat' may help adapting
     data to this requirement). The datasets must be global, i.e. spanning the
     whole Earth.
    d. In case variable and dimension names do not comply to these requirements
       you may use the CDO command (for variables): 
       'chname,namein,nameout filein.nc fileout.nc' 
       and the NCO command (for dimensions): 
       'ncrename -h -d namein,nameout -v .namein,nameout filename.nc'
    e. Create a folder in the input directory for each model you need to
       analyse, named after the model.

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

OUTPUT

The output directory contains the following NetCDF files:
    - (output directory):
        atmos_transp_mean_<model_name>.nc
        latent_transp_mean_<model_name>.nc
        ocean_transp_mean_<model_name>.nc
        total_transp_mean_<model_name>.nc
        wmb_transp_mean_<model_name>.nc

        contain annual mean meridional sections of heat transports in the
        atmosphere, oceans, and as a total; latent energy transports and water
        mass transports;

    - (output directory)/<model_name>:
        <model-name>_atmb.nc
        (<model-name>_latent.nc; if wat is set to true)
        <model-name>_surb.nc
        <model-name>_toab.nc
        (<model-name>_wmb.nc; is wat is set to true)

        contain annual mean 2D fields of energy budget, latent heat and water
        mass budgets;

        <model-name>_barocEff.nc

        contains the evolution of annual mean baroclinic efficiency
        (Lucarini et al., 2011).

        (if entr is set to true):
        <model-name>_evap_entr.nc (if met is set to 2 or 3)
        <model-name>_horizEntropy.nc (if met is set to 1 or 3)
        <model-name>_pot_drop_entr.nc (if met is set to 2 or 3)
        <model-name>_rain_entr.nc (if met is set to 2 or 3)
        <model-name>_sens_entr.nc (if met is set to 2 or 3)
        <model-name>_snow_entr.nc (if met is set to 2 or 3)
        <model-name>_snowmelt_entr.nc (if met is set to 2 or 3)
        <model-name>_verticalEntropy.nc (if met is set to 1 or 3)
        contain the evolution of annual mean components of the material entropy
        production.

    - (plots directory):
        meridional_transp.png: contains the model inter-comparison of total,
        atmospheric and oceanic of meridional sections in zonally averaged
        meridional heat transports;
        scatters_summary.png: contains the scatter plots of
        model intercomparisons of various metrics retrieved in the program;
        scattes_variability: contains scatter plots of model intercomparisons
        between TOA, atmospheric and surface global mean energy budgets and
        their inter-annual variability;

    - (plots directory)/<model-name>:
        <model-name>_atmb_timeser.png: the atmospheric budget annual mean
        global and hemispheric time series;
        <model-name>_energy_climap.png: the TOA, atmospheric and surface
        climatological mean fields;
        <model-name>_heat_transp.png: the meridional sections of total,
        atmospheric and oceanic meridional heat transports (implied from energy
        budgets);
        <model-name>_latent_climap.png: the climatological mean latent heat
        field;
        <model-name>_latent_timeser.png: the latent heat annual mean global and
        hemispheric evolutions;
        <model-name>_latent_transp.png: the meridional section of annual mean
        meridional latent heat transport;
        <model-name>_scatpeak.png: the scatter plots of atmospheric vs. oceanic
        peak magnitude in both hemispheres;
        <model-name>_sevap_climap.png: the annual mean field of material
        entropy production due to evaporation;
        <model-name>_smelt_climap.png: the annual mean field of material
        entropy production due to snow melting;
        <model-name>_spotp_climap.png: the annual mean field of material
        entropy production due to potential energy of the droplet;
        <model-name>_srain_climap.png: the annual mean field of material
        entropy production due to rainfall precipitation;
        <model-name>_ssens_climap.png: the annual mean field of material
        entropy production due to sensible heat fluxes;
        <model-name>_ssnow_climap.png: the annual mean field of material
        entropy production due to snowfall precipitation;
        <model-name>_surb_timeser.png: the surface budget annual mean
        global and hemispheric time series;
        <model-name>_sver_climap.png: the annual mean field of vertical
        material entropy production through the indirect method;
        <model-name>_toab_timeser.png: the TOA budget annual mean
        global and hemispheric time series;
        <model-name>_wmb_climap.png: the climatological mean water mass budget
        field;
        <model-name>_wmb_timeser.png: the water mass annual mean global and
        hemispheric evolutions;
        <model-name>_wmb_transp.png: the meridional section of annual mean
        meridional water mass transport;

    - (plots directory)/<model-name>/LEC_results:
        <model-name>_<year>_lec_diagram.png: the flux diagram for the annual
        mean LEC cycle in a specific year;
        <model-name>_<year>_lec_table.txt: the table containing the storage and
        conversion terms for the annual mean LEC cycle in a specific year;

#############################################################################

2019-04-29: a stand-alone version of TheDiaTo v1.0 is branched from ESMValTool
            v2.0b repository;

#############################################################################
"""

import logging
import os
import warnings
import glob
import numpy as np

# Locally used modules
from namelist import direc, models, flagin, logfile
import computations, lorenz_cycle, mkthe, plot_script

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
logging.basicConfig(filename=logfile, level=logging.INFO)
logger = logging.getLogger(__file__)

lorenz = lorenz_cycle
comp = computations
logger.info('Entering the diagnostic tool')
# Load paths
logger.info('Work directory: %s \n', direc[2])
logger.info('Plot directory: %s \n', direc[1])
plotsmod = plot_script
logger.info('model_names')
flags = [flagin[1], flagin[3], flagin[4]]
# Initialize multi-model arrays
modnum = len(models)
te_all = np.zeros(modnum)
toab_all = np.zeros([modnum, 2])
toab_oc_all = np.zeros(modnum)
toab_la_all = np.zeros(modnum)
atmb_all = np.zeros([modnum, 2])
atmb_oc_all = np.zeros(modnum)
atmb_la_all = np.zeros(modnum)
surb_all = np.zeros([modnum, 2])
surb_oc_all = np.zeros(modnum)
surb_la_all = np.zeros(modnum)
wmb_all = np.zeros([modnum, 2])
wmb_oc_all = np.zeros(modnum)
wmb_la_all = np.zeros(modnum)
latent_all = np.zeros([modnum, 2])
latent_oc_all = np.zeros(modnum)
latent_la_all = np.zeros(modnum)
baroc_eff_all = np.zeros(modnum)
lec_all = np.zeros([modnum, 2])
horzentr_all = np.zeros([modnum, 2])
vertentr_all = np.zeros([modnum, 2])
matentr_all = np.zeros([modnum, 2])
irrevers_all = np.zeros(modnum)
diffentr_all = np.zeros([modnum, 2])
logger.info("Entering main loop\n")
i_m = 0
for model in models:
    # Load paths to individual models output and plotting directories
    idir = os.path.join(direc[0], model)
    wdir = os.path.join(direc[2], model)
    pdir = os.path.join(direc[1], model)
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    if not os.path.exists(pdir):
        os.makedirs(pdir)
    # Reading file names for the specific model
    logger.info('Processing model: %s \n', model)
    filenames = [f for f in glob.glob(idir + "/*.nc", recursive=True)]
    filenames.sort()
    rlds_file = filenames[6]
    rlus_file = filenames[7]
    rsds_file = filenames[9]
    rsus_file = filenames[11]
    ta_file = filenames[13]
    ts_file = filenames[15]
    aux_file = wdir + '/aux.nc'
    te_ymm_file, te_gmean_constant, _, _ = mkthe.init_mkthe(
        model, wdir, filenames, flags)
    te_all[i_m] = te_gmean_constant
    logger.info('Computing energy budgets\n')
    eb_gmean, eb_file, toab_ymm_file = comp.budgets(
        model, wdir, aux_file, filenames)
    toab_all[i_m, 0] = np.nanmean(eb_gmean[0])
    toab_all[i_m, 1] = np.nanstd(eb_gmean[0])
    atmb_all[i_m, 0] = np.nanmean(eb_gmean[1])
    atmb_all[i_m, 1] = np.nanstd(eb_gmean[1])
    surb_all[i_m, 0] = np.nanmean(eb_gmean[2])
    surb_all[i_m, 1] = np.nanstd(eb_gmean[2])
    logger.info('Global mean emission temperature: %s\n',
                te_gmean_constant)
    logger.info('TOA energy budget: %s\n', toab_all[i_m, 0])
    logger.info('Atmospheric energy budget: %s\n', atmb_all[i_m, 0])
    logger.info('Surface energy budget: %s\n', surb_all[i_m, 0])
    logger.info('Done\n')
    baroc_eff_all[i_m] = comp.baroceff(model, wdir, aux_file,
                                       toab_ymm_file, te_ymm_file)
    logger.info('Baroclinic efficiency (Lucarini et al., 2011): %s\n',
                baroc_eff_all[i_m])
    logger.info('Running the plotting module for the budgets\n')
    plotsmod.balances(
        direc[2], pdir, [eb_file[0], eb_file[1], eb_file[2]],
        ['toab', 'atmb', 'surb'], model)
    logger.info('Done\n')
    # Water mass budget
    if flagin[1] == 'True':
        logger.info('Computing water mass and latent energy budgets\n')
        _, _, _, aux_list = mkthe.init_mkthe(model, wdir, filenames, flags)
        wm_gmean, wm_file = comp.wmbudg(model, wdir, aux_file, filenames,
                                        aux_list)
        wmb_all[i_m, 0] = np.nanmean(wm_gmean[0])
        wmb_all[i_m, 1] = np.nanstd(wm_gmean[0])
        logger.info('Water mass budget: %s\n', wmb_all[i_m, 0])
        latent_all[i_m, 0] = np.nanmean(wm_gmean[1])
        latent_all[i_m, 1] = np.nanstd(wm_gmean[1])
        logger.info('Latent energy budget: %s\n', latent_all[i_m, 0])
        logger.info('Done\n')
        logger.info('Plotting the water mass and latent energy budgets\n')
        plotsmod.balances(direc[2], pdir, [wm_file[0], wm_file[1]],
                          ['wmb', 'latent'], model)
        logger.info('Done\n')
        for filen in aux_list:
            os.remove(filen)
    if flagin[0] == 'True':
        lsmdir = os.path.join(direc[0], 'ls_mask')
        sftlf_fx = lsmdir + '/lsmask_{}.nc'.format(model)
        logger.info('Computing energy budgets over land and oceans\n')
        toab_oc_gmean, toab_la_gmean = comp.landoc_budg(
            model, wdir, eb_file[0], sftlf_fx, 'toab')
        toab_oc_all[i_m] = toab_oc_gmean
        toab_la_all[i_m] = toab_la_gmean
        logger.info('TOA energy budget over oceans: %s\n', toab_oc_gmean)
        logger.info('TOA energy budget over land: %s\n', toab_la_gmean)
        atmb_oc_gmean, atmb_la_gmean = comp.landoc_budg(
            model, wdir, eb_file[1], sftlf_fx, 'atmb')
        atmb_oc_all[i_m] = atmb_oc_gmean
        atmb_la_all[i_m] = atmb_la_gmean
        logger.info('Atmospheric energy budget over oceans: %s\n',
                    atmb_oc_gmean)
        logger.info('Atmospheric energy budget over land: %s\n',
                    atmb_la_gmean)
        surb_oc_gmean, surb_la_gmean = comp.landoc_budg(
            model, wdir, eb_file[2], sftlf_fx, 'surb')
        surb_oc_all[i_m] = surb_oc_gmean
        surb_la_all[i_m] = surb_la_gmean
        logger.info('Surface energy budget over oceans: %s\n',
                    surb_oc_gmean)
        logger.info('Surface energy budget over land: %s\n', surb_la_gmean)
        logger.info('Done\n')
        if flagin[1] == 'True':
            logger.info('Computing water mass and latent energy'
                        ' budgets over land and oceans\n')
            wmb_oc_gmean, wmb_la_gmean = comp.landoc_budg(
                model, wdir, wm_file[0], sftlf_fx, 'wmb')
            wmb_oc_all[i_m] = wmb_oc_gmean
            wmb_la_all[i_m] = wmb_la_gmean
            logger.info('Water mass budget over oceans: %s\n',
                        wmb_oc_gmean)
            logger.info('Water mass budget over land: %s\n', wmb_la_gmean)
            latent_oc_gmean, latent_la_gmean = comp.landoc_budg(
                model, wdir, wm_file[1], sftlf_fx, 'latent')
            latent_oc_all[i_m] = latent_oc_gmean
            latent_la_all[i_m] = latent_la_gmean
            logger.info('Latent energy budget over oceans: %s\n',
                        latent_oc_gmean)
            logger.info('Latent energy budget over land: %s\n',
                        latent_la_gmean)
            logger.info('Done\n')
    if flagin[2] == 'True':
        logger.info('Computation of the Lorenz Energy '
                    'Cycle (year by year)\n')
        ldir = os.path.join(pdir, 'LEC_results')
        if not os.path.exists(ldir):
            os.makedirs(ldir)
        lect = lorenz.preproc_lec(model, wdir, ldir, filenames)
        lec_all[i_m, 0] = np.nanmean(lect)
        lec_all[i_m, 1] = np.nanstd(lect)
        logger.info(
            'Intensity of the annual mean Lorenz Energy '
            'Cycle: %s\n', lec_all[i_m, 0])
        logger.info('Done\n')
    else:
        lect = np.repeat(2.0, len(eb_gmean[0]))
        lec_all[i_m, 0] = 2.0
        lec_all[i_m, 1] = 0.2
    if flagin[3] == 'True':
        if flagin[4] in {'1', '3'}:
            _, _, te_file, _ = mkthe.init_mkthe(model, wdir, filenames,
                                                flags)
            logger.info('Computation of the material entropy production '
                        'with the indirect method\n')
            indentr_list = [
                rlds_file, rlus_file, rsds_file, rsus_file, te_file,
                eb_file[0], ts_file
            ]
            horz_mn, vert_mn, horzentr_file, vertentr_file = comp.indentr(
                model, wdir, indentr_list, aux_file, eb_gmean[0])
            horzentr_all[i_m, 0] = np.nanmean(horz_mn)
            horzentr_all[i_m, 1] = np.nanstd(horz_mn)
            vertentr_all[i_m, 0] = np.nanmean(vert_mn)
            vertentr_all[i_m, 1] = np.nanstd(vert_mn)
            logger.info(
                'Horizontal component of the material entropy '
                'production: %s\n', horzentr_all[i_m, 0])
            logger.info(
                'Vertical component of the material entropy '
                'production: %s\n', vertentr_all[i_m, 0])
            logger.info('Done\n')
            logger.info('Running the plotting module for the material '
                        'entropy production (indirect method)\n')
            plotsmod.entropy(pdir, vertentr_file, 'sver',
                             'Vertical entropy production', model)
            os.remove(te_file)
            logger.info('Done\n')
        if flagin[4] in {'2', '3'}:
            matentr, irrevers, entr_list = comp.direntr(
                logger, model, wdir, filenames, aux_file, lect,
                flagin[2], flags)
            matentr_all[i_m, 0] = matentr
            if flagin[4] in {'3'}:
                diffentr = (float(np.nanmean(vert_mn)) + float(
                    np.nanmean(horz_mn)) - matentr)
                logger.info('Difference between the two '
                            'methods: %s\n', diffentr)
                diffentr_all[i_m, 0] = diffentr
            logger.info('Degree of irreversibility of the '
                        'system: %s\n', irrevers)
            irrevers_all[i_m] = irrevers
            logger.info('Running the plotting module for the material '
                        'entropy production (direct method)\n')
            plotsmod.init_plotentr(model, pdir, entr_list)
            logger.info('Done\n')
    os.remove(te_ymm_file)
    logger.info('Done for model: %s \n', model)
    i_m = i_m + 1
logger.info('I will now start multi-model plots')
logger.info('Meridional heat transports\n')
plotsmod.plot_mm_transp(models, direc[2], direc[1])
logger.info('Scatter plots')
summary_varlist = [
    atmb_all, baroc_eff_all, horzentr_all, lec_all, matentr_all, te_all,
    toab_all, vertentr_all
]
plotsmod.plot_mm_summaryscat(direc[1], summary_varlist)
logger.info('Scatter plots for inter-annual variability of'
            ' some quantities')
eb_list = [toab_all, atmb_all, surb_all]
plotsmod.plot_mm_ebscatter(direc[1], eb_list)
logger.info("The diagnostic has finished. Now closing...\n")
