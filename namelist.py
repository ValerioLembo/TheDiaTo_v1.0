"""User options for TheDiaTo v1.0.

Write here the path to the input, output and work directories, the model names
that are analysed. Switch on and off here the different modules. Set the name
for the output log containing the values of the global metrics.

############################################################################
2019-04-29: the script has been created.
############################################################################

Author: Valerio Lembo, Meteorological Institute, University of Hamburg
E-mail: valerio.lembo@uni-hamburg.de .
"""

import datetime

now = datetime.datetime.now()

# Write here the path to the input directory, where the model outputs are
# contained, the plots directory, where the plots and the tables of the LEC are
# saved, the working directory, where the output NetCDF datasets are stored.
#
# The model outputs for each model shall be stored in a different subfolder in
# the input parent directory. This subfolder must be named after the model.
#
# A subfolder with the name of the model is automatically created in the plots
# and working directories.
idir_up = '/work/um0005/u234097/ESMV/modeldata_standalone'
pdir_up = '/work/um0005/u234097/ESMV/plots_standalone'
wdir_up = '/work/um0005/u234097/ESMV/output_standalone'
direc = [idir_up, pdir_up, wdir_up]

# Write the model names you want to analyse in this list, as strings and
# separated by commas.
models = ['MPI-ESM-P']

# Set the flags for the modules to 'True' or 'False'. Set the 'met' flag to
# '1' for computing the MEP with the indirect method, '2' for the direct method
# , '3' for both.
lsm = 'False' # Flag for the land-ocean computations
wat = 'False' # Flag for the water and latent energy budget
lec = 'True' # Flag for the LEC
entr = 'True' # Flag for the MEP
met = '3' # Option for the MEP method

flagin = [lsm, wat, lec, entr, met]
date = now.isoformat()
logfile = 'log_{}.txt'.format(date) # Put the desired name for the log here