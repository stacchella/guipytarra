"""
Helpher functions to run guitarra.

    convert_str:
        Function converts some input variable to a Fortran-like
        string.

    setup_input_files:
        Function converts input parameters to two input files
        that are read by guitarra.

    run_guitarra:
        Main function that runs guitarra.

"""

# --------------
# import modules
# --------------

import os


# --------------
# define functions
# --------------

def convert_str(input_obj):
    """
    Function converts some input variable to a Fortran-like
    string.

    Parameters
    ----------
    input_obj : any type
        This variable is converted to a string that can be
        read by Fotran.

    """
    if isinstance(input_obj, str):
        return(input_obj)
    elif input_obj is True:
        return('1')
    elif (input_obj is None):
        return('none')
    elif input_obj is False:
        return('0')
    else:
        return(str(input_obj))


def setup_input_files(parameter_dictionary):
    """
    Function converts input parameters to two input files
    that are read by guitarra.

    Parameters
    ----------
    parameter_dictionary : dict
        Dictionary that contains all parameters needed to run
        guitarra.

    """
    # parameter file that is conveyed in batch mode
    list_parameters_1 = ['filename_param2', 'idither', 'filename', 'noise_name', 'ra0', 'dec0', 'new_ra',
                         'new_dec', 'dx(idither)', 'dy(idither)', 'sca_id', 'indx', 'icat_f',
                         'star_catalogue', 'galaxy_catalogue', 'filters_in_cat']
    file_1 = open(parameter_dictionary['filename_param1'], 'w')
    for ii_key in list_parameters_1:
        file_1.write(convert_str(parameter_dictionary[ii_key]) + '\n')
    file_1.close()
    # parameter file that is conveyed in the background
    list_parameters_2 = ['verbose', 'brain_dead_test', 'module', 'mode', 'ngroups', 'subarray',
                  'colcornr', 'rowcornr', 'naxis1', 'naxis2', 'number_primary', 'number_subpixel',
                  'camera', 'primary_dither', 'subpixel_dither', 'ra0', 'dec0', 'pa_degrees',
                  'include_ktc', 'include_dark', 'include_readnoise', 'include_reference', 'include_non_linear', 'include_latents',
                  'include_1_over_f', 'include_cr', 'cr_mode', 'include_bg', 'bkg_mode', 'zodiacal_scale_factor',
                  'include_stars', 'nstars', 'star_catalogue', 'include_galaxies', 'ngal', 'galaxy_catalogue',
                  'nf', 'use_filter', 'cat_filter']
    file_2 = open(parameter_dictionary['filename_param2'], 'w')
    for ii_key in list_parameters_2:
        file_2.write(convert_str(parameter_dictionary[ii_key]) + '\n')
    file_2.close()


def run_guitarra(parameter_dictionary):
    """
    Main function that runs guitarra.

    Parameters
    ----------
    parameter_dictionary : dict
        Dictionary that contains all parameters needed to run
        guitarra.

    """
    # setup input files
    setup_input_files(parameter_dictionary)
    # run guitarra
    command = os.environ['GUITARRA_HOME']+'/src/fortran/guitarra < ' + parameter_dictionary['filename_param1']
    os.system(command)
