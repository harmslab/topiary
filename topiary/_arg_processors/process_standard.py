__description__ = \
"""
Functions for processing/checking standard variable inputs (bool, float, int,
etc.)
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-27"

import numpy as np

def process_bool(value,variable_name=None):
    """
    Process a bool argument and do error checking.

    Parameters
    ----------
        value: input value to check/process
        variable_name: name of variable (string, for error message)

    Return
    ------
        validated value
    """

    try:

        # See if this is an iterable
        if hasattr(value,"__iter__"):
            raise ValueError

        # See if this is a naked type
        if type(value) is type:
            raise ValueError

        if value != 0:
            if not np.isclose(round(value,0)/value,1):
                raise ValueError

        value = bool(int(value))

    except (TypeError,ValueError):

        if variable_name is not None:
            err = f"\n{variable_name} '{value}' must be True or False.\n\n"
        else:
            err = f"\n'{value}' must be True or False.\n\n"
        err += "\n\n"

        raise ValueError(value)

    return value


def process_float(value,
                  variable_name=None,
                  minimum_allowed=-np.inf,
                  maximum_allowed=np.inf,
                  minimum_inclusive=True,
                  maximum_inclusive=True):
    """
    Process a float argument and do error checking.

    Parameters
    ----------
        value: input value to check/process
        variable_name: name of variable (string, for error message)
        minimum_allowed: minimum allowable value for the variable
        maximum_allowed: maximum allowable value for the variable
        minimum_inclusive: whether lower bound is inclusive
        maximum_inclusive: whether upper bound is inclusive

    Return
    ------
        validated value
    """

    try:

        # Try to cast as string to an integer
        if type(value) is str:
            value = float(value)

        # See if this is an iterable
        if hasattr(value,"__iter__"):
            raise ValueError

        # See if this is a naked type
        if type(value) is type:
            raise ValueError

        value = float(value)

        if np.isnan(value):
            raise ValueError

        if minimum_inclusive:
            if value < minimum_allowed:
                raise ValueError
        else:
            if value <= minimum_allowed:
                raise ValueError

        if maximum_inclusive:
            if value > maximum_allowed:
                raise ValueError
        else:
            if value >= maximum_allowed:
                raise ValueError

    except (ValueError,TypeError):

        if minimum_inclusive:
            min_o = "<="
        else:
            min_o = "<"

        if maximum_inclusive:
            max_o = "<="
        else:
            max_o = "<"

        bounds = f"{minimum_allowed} {min_o} {variable_name} {max_o} {maximum_allowed}"

        if variable_name is not None:
            err = f"\n{variable_name} '{value}' must be a float:\n\n"
        else:
            err = f"\n'{value}' must be a float:\n\n"

        if not (minimum_allowed is None and maximum_allowed is None):
            err += bounds

        err += "\n\n"

        raise ValueError(err)

    return value

def process_int(value,
                variable_name=None,
                minimum_allowed=None,
                maximum_allowed=None,
                minimum_inclusive=True,
                maximum_inclusive=True):
    """
    Process an integer argument and do error checking.

    Parameters
    ----------
        value: input value to check/process
        variable_name: name of variable (string, for error message)
        minimum_allowed: minimum allowable value for the variable
        maximum_allowed: maximum allowable value for the variable
        minimum_inclusive: whether lower bound is inclusive
        maximum_inclusive: whether upper bound is inclusive

    Return
    ------
        validated value
    """


    try:

        # Try to cast as string to an integer
        if type(value) is str:
            value = int(value)

        # See if this is an iterable
        if hasattr(value,"__iter__"):
            raise ValueError

        # See if this is a naked type
        if type(value) is type:
            raise ValueError

        # If this is a float to int cast, make sure it does not have decimal
        floored = np.floor(value)
        if value != floored:
            raise ValueError

        # Make int cast
        value = int(value)

        if minimum_allowed is not None:
            if minimum_inclusive:
                if value < minimum_allowed:
                    raise ValueError
            else:
                if value <= minimum_allowed:
                    raise ValueError

        if maximum_allowed is not None:
            if maximum_inclusive:
                if value > maximum_allowed:
                    raise ValueError
            else:
                if value >= maximum_allowed:
                    raise ValueError

    except (ValueError,TypeError,OverflowError):

        if minimum_inclusive:
            min_o = "<="
        else:
            min_o = "<"

        if maximum_inclusive:
            max_o = "<="
        else:
            max_o = "<"

        bounds = f"{minimum_allowed} {min_o} {variable_name} {max_o} {maximum_allowed}"

        if variable_name is not None:
            err = f"\n{variable_name} '{value}' must be an integer:\n\n"
        else:
            err = f"\n'{value}' must be an integer:\n\n"

        if not (minimum_allowed is None and maximum_allowed is None):
            err += bounds

        err += "\n\n"

        raise ValueError(err)

    return value

def process_iter(value,
                 variable_name=None,
                 required_iter_type=None,
                 required_value_type=None,
                 minimum_allowed=None,
                 maximum_allowed=None,
                 minimum_inclusive=True,
                 maximum_inclusive=True,
                 is_not_type=None):
    """
    Process an iteratble and do error checking.

    NOTE: This function does *not* check required_value_type for numpy arrays.

    Parameters
    ----------
        value: input value to check/process
        variable_name: name of variable (string, for error message)
        required_iter_type: if not None, validate that the iterable has the
                            specified type
        required_value_type: if not None, validate that every value in the
                             iterable has the specified type
        minimum_allowed: minimum allowable length for the iterable
        maximum_allowed: maximum allowable length for the iterable
        minimum_inclusive: whether lower bound is inclusive
        maximum_inclusive: whether upper bound is inclusive
        is_not_type: type (or list of types) to reject

    Return
    ------
        validated value
    """

    # Set up error message
    if variable_name is not None:
        err_base = f"\nvariable_name '{value}' "
    else:
        err_base = f"\n'{value}' "

    # Make sure it's an iterable
    if not hasattr(value,"__iter__"):
        err = err_base + "must be list-like\n"
        raise ValueError(err)

    # See if this is a naked type
    if type(value) is type:
        err = err_base + "must not be a type\n"
        raise ValueError(err)

    # If requested, make sure it's the required type
    if required_iter_type is not None:
        if type(value) is not required_iter_type:
            err = err_base + f"must be type {required_iter_type}\n"
            raise ValueError(err)

    # If reuqested, make sure the values are of the right type
    if required_value_type is not None:

        # Skip check for iterable with no values
        if len(value) > 0:

            # Only do check iterables that are not numpy arrays
            if type(value) is not np.ndarray:
                types = list(set([type(v) for v in value]))
                if len(types) != 1 or types[0] is not required_value_type:
                    err = err_base + f"all entries must have type {required_value_type}\n"
                    raise ValueError(err)

    # If requested, make sure the iterable has appropriate dimensions
    try:
        if minimum_allowed is not None:
            if minimum_inclusive:
                if len(value) < minimum_allowed:
                    raise ValueError
            else:
                if len(value) <= minimum_allowed:
                    raise ValueError

        if maximum_allowed is not None:
            if maximum_inclusive:
                if len(value) > maximum_allowed:
                    raise ValueError
            else:
                if len(value) >= maximum_allowed:
                    raise ValueError
    except ValueError:

        if minimum_inclusive:
            min_o = "<="
        else:
            min_o = "<"

        if maximum_inclusive:
            max_o = "<="
        else:
            max_o = "<"

        bounds = f"{minimum_allowed} {min_o} {variable_name} {max_o} {maximum_allowed}"

        err = err_base + f"must have length:\n\n{bounds}\n\n"
        raise ValueError(err)

    # Make sure iterable does not have a specified excluded type
    if is_not_type is not None:
        if type(is_not_type) is type:
            is_not_type = [is_not_type]

        if type(value) in is_not_type:
            bad_types = ",".join([f"{v}" for v in is_not_type])

            err = err_base + f"must not be of type {bad_types}\n\n"
            raise ValueError(err)

    return value
