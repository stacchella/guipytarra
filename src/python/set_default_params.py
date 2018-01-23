'''
This module sets the default parameters for guipytarra.
'''

# import modules

import os
import pickle


# set up parameter dictionary

run_params = {'idither': 90,  # something to do with dither position? (f.e. 90)
              'filename': 'sim_cube_F200W_487_090.fits',  # output file name (f.e. )
              'filename_param1': 'batch_input.input',  # filename of input via batch
              'filename_param2': 'simulator.params',  # filename of input via background
              'noise_name': 'TBD',  # something to do with noise? (f.e. None)  => to do with add_one_over_f_noise
              # image RA and DEC information
#twice              'ra0': 53.152793,  # some initial ra?
#twice              'dec0': -27.769994,    # some initial dec?
              'new_ra': 53.152793,  # some new_ra
              'new_dec': -27.769994,  # some new dec?
              'dx(idither)': 53.152793,  # ?
              'dy(idither)': -27.769994,  # ?
              # detector, filter and catalog input
              'sca_id': 487,  # ID of SCA (=> but already written by read_parameters => modules?) (f.e. 487)
              'indx': 10,  # related to filter in 'nircam_calib.list' and 'psf.list'  (f.e. 10 for F200W)
              'icat_f': 8,  # something filter (f.e. 8 for ?) in catalogue?
#twice              'star_catalogue': 'none_487_090.cat',  # catalog of stars, see below (f.e. none_487_090.cat)
#twice              'galaxy_catalogue': 'mock_2017_11_03_487_090.cat',  # catalog of galaxies, see below (f.e. mock_2017_11_03_487_090.cat)
              #
              'verbose': 1,  # verbose (f.e. 1)
              'brain_dead_test': False,  # brain_dead_test (1=yes)
              'module': 'NRCALL',  # module (f.e. NRCALL)  show other options
              'mode': 'deep8',  # exposure mode (f.e. deep8)
              'ngroups': 7,  # number of groups per ramp  (f.e. 7)
              'subarray': 'FULL',  # Subarray mode (f.e. FULL)
              'colcornr': 0,  # lower right corner X coordinate (f.e. 0)
              'rowcornr': 0,  # Lower right corner Y coordinate (f.e. 0)
              'naxis1': 2048,  # number of X pixels (f.e. 2048)
              'naxis2': 2048,  # number of Y pixels (f.e. 2048)
              'number_primary': 1,  # number of primary dithers (f.e. 1)
              'number_subpixel': 9,  # number of subpixel dithers (f.e. 9)
              'camera': 'short',  # camera (f.e. short)
              'primary_dither': 'intrasca',  # primary dither (f.e. intrasca)
              'subpixel_dither': 'small',  # subpixel dither  (f.e. small)
              'ra0': 53.1777830,  # RA (degrees) of center (f.e. 53.1777830)  ????
              'dec0': -27.7999000,  # Dec (degrees) of center (f.e. -27.7999000) ????
              'pa_degrees': 295.0264750,  # Position angle (degrees) of NIRCam short axis, N to E (f.e. 295.0264750)
              'include_ktc': True,  # include ktc
              'include_dark': True,  # include dark
              'include_readnoise': True,  # include readnoise
              'include_reference': True,  # include reference
              'include_non_linear': True,  # include non-linearity
              'include_latents': False,  # include latents
              'include_1_over_f': False,  # include 1/f noise
              'include_cr': False,  # include cosmic rays
              'cr_mode': 1,  # cosmic ray mode: (0) ACS; (1) quiet Sun (2); active Sun; (3) flare
              'include_bg': True,  # include zodiacal background
              'bkg_mode': 3,  # mode (0)  MJR; (2) STScI model (3) SED convolution
              'zodiacal_scale_factor': 1,  # verbose (f.e. 1)
              'include_stars': False,  # include stars
              'nstars': 0,  # number of stars (f.e. 0)
              'star_catalogue': None,  # galaxy catalogue (f.e. none)
              'include_galaxies': True,  # include galaxies
              'ngal': 10,  # number of galaxies (f.e. 10)
              'galaxy_catalogue': 'mock_2017_11_03_487_090.cat',  # galaxy catalogue (f.e. mock_2017_11_03.cat)
              'filters_in_cat': 20,   # number of filters in catalog
              'nf': 1,  # number of filters to use (f.e. 2)
              'use_filter': 10,  # filter name (f.e. 10 => F200W, 18 => F356W)
              'cat_filter': 8,  # filter number in catalogue
               }

# save parameter dictionary

with open(os.environ['GUITARRA_HOME'] + 'data/default_params.pkl', 'wb') as f:
    pickle.dump(run_params, f, pickle.HIGHEST_PROTOCOL)


