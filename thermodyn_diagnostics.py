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

For additional information and user manual see README.md

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

list_basic=[
    'hfls_','hfss_','rlds_','rlus_','rlut_','rsds_','rsdt_','rsus_','rsut_']
list_wat=['pr_','prsn_']
list_lec=['ta_','tas_','ua_','uas_','va_','vas_','wap_']
list_indentr=['ts_']
list_direntr=['hus_','pr_','prsn_','ps_','ts_']
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
flags = [flagin[1], flagin[2], flagin[3], flagin[4]]
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
    for i in list_basic:
        for name in filenames:
            if i in name:
                exec("%sfile = '%s'" % (i,name))
    #rlds_file = filenames[6]
    #rlus_file = filenames[7]
    #rsds_file = filenames[9]
    #rsus_file = filenames[11]
    #ts_file = filenames[15]
    aux_file = wdir + '/aux.nc'
    te_ymm_file, te_gmean_constant, _, _ = mkthe.init_mkthe(
        model, wdir, rlut_file, flags)
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
            plotsmod.entropy(pdir, horzentr_file, 'shor',
                             'Horizontal entropy production', model)
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
