'''
This module sets the default parameters for guipytarra.
'''

# import modules

import os
import pickle


# set up parameter dictionary

run_params = {'verbose': 1,                # verbose (f.e. 1)
              'brain_dead_test': False,    # brain_dead_test (1=yes)
              # define names of files
              'idither': 1,  # number of dither position
              'filename': 'sim_cube_filter_sca_idither.fits',  # output file name
              'filename_param': 'batch_input.input',           # filename of input via batch
              'noise_name': 'noise_file.txt',         # noise map of one over f noise
              'zodifile': 'bkg_file.txt',  # name of file containing background SED for observation date
              # image RA and DEC information
              'ra0': 53.17778300,          # RA of pointing
              'dec0': -27.79990000,        # DEC of pointing
              # detector, filter and catalog input
              'sca_id': 490,               # ID of SCA to use  (f.e. 487)
              'filter_index': 14,          # index of filter in 'nircam_calib.list' and 'psf.list'  (f.e. 10 for F200W)
              'include_stars': False,      # include stars
              'star_catalogue': 'star_cat_filename.cat',    # catalog of stars, see below (f.e. none_487_090.cat)
              'include_galaxies': True,    # include galaxies
              'galaxy_catalogue': 'galaxy_cat_filename.cat',  # catalog of galaxies, see below (f.e. mock_2017_11_03_487_090.cat)
              'include_cloned_galaxies': False,  # include cloned galaxies
              'filter_in_cat': 8,          # number of filters contained in source catalogues
              'icat_f': 5,                 # index of filter to use
              'psf_add': True,             # convolve with PSF
              # aperture
              'apername': 'BLONG',         # aperture name
              'readpatt': 'deep8',         # read out pattern
              'ngroups': 7,                # number of groups per ramp  (f.e. 7)
              'subarray': 'FULL',          # subarray mode (f.e. FULL)
              'substrt1': 0,               # lower right corner X coordinate (f.e. 0)
              'substrt2': 0,               # lower right corner Y coordinate (f.e. 0)
              'subsize1': 2048,            # number of X pixels (f.e. 2048)
              'subsize2': 2048,            # number of Y pixels (f.e. 2048)
              'pa_degrees': 295.026475,    # Position angle (degrees) of NIRCam short axis, N to E (f.e. 295.0264750)
              # noise to include
              'noiseless': False,          # run without adding noise
              'ipc_add': True,             # add interpixel capacitance
              'include_ktc': True,         # include ktc
              'include_dark': True,        # include dark
              'include_readnoise': True,   # include readnoise
              'include_reference': True,   # include reference
              'include_non_linear': True,  # include non-linearity
              'include_latents': False,    # include latents
              'include_1_over_f': False,   # include 1/f noise
              'include_cr': False,         # include cosmic rays
              'cr_mode': 1,                # cosmic ray mode: (0) ACS; (1) quiet Sun (2); active Sun; (3) flare
              'include_bg': True,          # include zodiacal background
              }

# save parameter dictionary

with open(os.environ['GUITARRA_HOME'] + 'data/default_params.pkl', 'wb') as f:
    pickle.dump(run_params, f, pickle.HIGHEST_PROTOCOL)


