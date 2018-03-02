"""
Helpher functions to run guitarra.

    convert_str:
        Function converts some input variable to a Fortran-like
        string.

    setup_input_file:
        Function converts input parameters to input file
        that is read by guitarra.

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


def setup_input_file(parameter_dictionary):
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
    list_parameters = ['verbose', 'brain_dead_test', 'idither', 'filename', 'noise_name', 'ra0', 'dec0', 'new_ra', 'new_dec',
                       'dx', 'dy', 'sca_id', 'indx', 'include_stars', 'star_catalogue', 'nstars', 'include_galaxies',
                       'galaxy_catalogue', 'ngal', 'include_cloned_galaxies', 'filter_in_cat', 'icat_f', 'zodifile',
                       'apername', 'readpatt', 'ngroups', 'subarray', 'substrt1', 'substrt2', 'subsize1', 'subsize2',
                       'pa_degrees', 'include_ktc', 'include_dark', 'include_readnoise', 'include_reference',
                       'include_non_linear', 'include_latents', 'include_1_over_f', 'include_cr', 'cr_mode', 'include_bg']
    file_batch = open(parameter_dictionary['filename_param'], 'w')
    for ii_key in list_parameters:
        file_batch.write(convert_str(parameter_dictionary[ii_key]) + '\n')
    file_batch.close()



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
    setup_input_file(parameter_dictionary)
    # run guitarra
    command = os.environ['GUITARRA_HOME']+'/src/fortran/guitarra < ' + parameter_dictionary['filename_param']
    os.system(command)
