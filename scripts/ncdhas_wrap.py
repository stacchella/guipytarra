'''
This script collects the input for guipytarra and
also runs it.
'''

# --------------
# import modules
# --------------

import os
from astropy.io import fits

print os.environ['NCDHAS_PATH']


# -------------------------------
# define file list name to reduce
# -------------------------------

file_name_list = ['/Users/sandrotacchella/ASTRO/JWST/img_simulator/test_arizona/sandro/sim_cube_F200W_481_001.fits']


# ------------
# define paths
# ------------

# Sandro Laptop
path_calibration_files = '../../../../../../Volumes/Tacchella/Work/Postdoc/JWST_GTO/cal/'

# Cluster
path_calibration_files = 'TBD'


# ----------------------------------
# define reduction specific features
# ----------------------------------

ncdhas_flags = ' +vb +ow +wi +ws -mr 2048 -zi +cbp +cs +cbs +cd +cl -rhf -rss -ipc -df 0 +cfg isimcv3 '

cfg_config = {}

cfg_config['481'] = path_calibration_files + 'SCAConfig/NRCA1_17004_SW_ISIMCV3.cfg'
cfg_config['482'] = path_calibration_files + 'SCAConfig/NRCA2_17006_SW_ISIMCV3.cfg'
cfg_config['483'] = path_calibration_files + 'SCAConfig/NRCA3_17012_SW_ISIMCV3.cfg'
cfg_config['484'] = path_calibration_files + 'SCAConfig/NRCA4_17048_SW_ISIMCV3.cfg'
cfg_config['485'] = path_calibration_files + 'SCAConfig/NRCA5_17158_LW_ISIMCV3.cfg'
cfg_config['486'] = path_calibration_files + 'SCAConfig/NRCB1_16991_SW_ISIMCV3.cfg'
cfg_config['487'] = path_calibration_files + 'SCAConfig/NRCB2_17005_SW_ISIMCV3.cfg'
cfg_config['488'] = path_calibration_files + 'SCAConfig/NRCB3_17011_SW_ISIMCV3.cfg'
cfg_config['489'] = path_calibration_files + 'SCAConfig/NRCB4_17047_SW_ISIMCV3.cfg'
cfg_config['490'] = path_calibration_files + 'SCAConfig/NRCB5_17161_LW_ISIMCV3.cfg'


# ------------
# define paths
# ------------

for ii in range(len(file_name_list)):
    print 'progress in %: ', 100.0*ii/len(file_name_list)
    hdu_list = fits.open(file_name_list[ii])
    sca_id = str(hdu_list[0].header['SCA_ID'])
    command = 'ncdhas ' + file_name_list[ii] + ncdhas_flags + '+cp ' + path_calibration_files + ' -P ' + cfg_config[sca_id]
    print 'executing... ', command
    os.system(command)





