import os
print os.environ['NCDHAS_PATH']


path_calibration_files = '../../../../../../Volumes/Tacchella/Work/Postdoc/JWST_GTO/cal/'

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


sca_id = '481'
file_name = '/Users/sandrotacchella/ASTRO/JWST/img_simulator/test_arizona/sandro/sim_cube_F200W_481_001.fits'

ncdhas_flags = ' +vb +ow +wi +ws -mr 2048 -zi +cbp +cs +cbs +cd +cl -ipc -df 1 +cfg isimcv3 +cp ' + path_calibration_files + ' -P ' + cfg_config[sca_id]
# i.e. no flat field correction, no IPC deconvolution
# maybe change -df 0 to -df 1?

print 'ncdhas ' + file_name + ncdhas_flags


