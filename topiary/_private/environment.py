"""
Functions for interacting with the operating system environment.
"""

import os

def load_env_variable(env_name,
                      check_function=None,
                      check_function_kwargs={}):
    """
    Load an environment variable, optionally parsing the type. 

    Parameters
    ----------
    env_name : str
        environment variable name
    check_function : function
        function to use to parse/check value. (i.e. check.check_int)
    check_function_kwargs : dict
        keyword arguments to pass to the check function

    Returns
    -------
    value : 
        environment variable value, typed appropriately. If env_name is not 
        defined, return None. 
    """

    value = None
    if env_name in os.environ:

        value = os.environ[env_name].strip()
        if check_function:
            try:
                value = check_function(value,**check_function_kwargs)
            except ValueError as error:
                err = f"environment variable {env_name} is defined by is not\n"
                err += f"the right type. Value: ${value}."

                raise ValueError(err) from error
    
    return value