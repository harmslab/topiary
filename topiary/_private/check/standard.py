"""
Functions to check/process `bool`, `float`, `int`, and iterable arguments in
topiary functions.
"""

import numpy as np

import re

def check_bool(value,variable_name=None):
    """
    Process a `bool` argument and do error checking.

    Parameters
    ----------
    value :
        input value to check/process
    variable_name : str
        name of variable (string, for error message)

    Returns
    -------
    bool
        validated/coerced bool

    Raises
    ------
    ValueError
        If value cannot be interpreted as a bool
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


def check_float(value,
                  variable_name=None,
                  minimum_allowed=-np.inf,
                  maximum_allowed=np.inf,
                  minimum_inclusive=True,
                  maximum_inclusive=True):
    """
    Process a `float` argument and do error checking.

    Parameters
    ----------
    value :
        input value to check/process
    variable_name : str
        name of variable (string, for error message)
    minimum_allowed : float, default=-np.inf
        minimum allowable value for the variable
    maximum_allowed : float, default=np.inf
        maximum allowable value for the variable
    minimum_inclusive : bool, default=True
        whether lower bound is inclusive
    maximum_inclusive : bool, default=True
        whether upper bound is inclusive

    Returns
    -------
    float
        validated/coerced float

    Raises
    ------
    ValueError
        If value cannot be interpreted as a float
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

def check_int(value,
                variable_name=None,
                minimum_allowed=None,
                maximum_allowed=None,
                minimum_inclusive=True,
                maximum_inclusive=True):
    """
    Process an `int` argument and do error checking.

    Parameters
    ----------
    value :
        input value to check/process
    variable_name : str
        name of variable (string, for error message)
    minimum_allowed : float, optional
        minimum allowable value for the variable
    maximum_allowed : float, optional
        maximum allowable value for the variable
    minimum_inclusive : bool, default=True
        whether lower bound is inclusive
    maximum_inclusive : bool, default=True
        whether upper bound is inclusive

    Returns
    -------
    int
        validated/coerced integer

    Raises
    ------
    ValueError
        If value cannot be interpreted as a int
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

def check_iter(value,
                 variable_name=None,
                 required_iter_type=None,
                 required_value_type=None,
                 minimum_allowed=None,
                 maximum_allowed=None,
                 minimum_inclusive=True,
                 maximum_inclusive=True,
                 is_not_type=None):
    """
    Process an iterable argument and do error checking.

    Parameters
    ----------
    value :
        input value to check/process
    variable_name : str
        name of variable (string, for error message)
    required_iter_type : type, optional, default=None
        If not None, validate that the iterable has the specified type
    required_value_type : type, optional, default=None
        if not None, validate that every value in the iterable has the
        specified type
    minimum_allowed : int, optional, default=None
        minimum allowable length for the iterable
    maximum_allowed : int, optional, default=None
        maximum allowable length for the iterable
    minimum_inclusive : bool, default=True
        whether lower bound is inclusive
    maximum_inclusive : bool, default=True
        whether upper bound is inclusive
    is_not_type : type, list, default=None
        type (or list of types) to reject

    Returns
    -------
    iterable
        validated/coerced iterable

    Notes
    -----
    This function does *not* check required_value_type for numpy arrays.

    Raises
    ------
    ValueError
        If value cannot be interpreted as an iterable of appropriate type
    """

    # Set up error message
    if variable_name is not None:
        err_base = f"\n{variable_name} = {value} "
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

def column_to_bool(column,column_name):
    """
    Convert a generic pandas column to bool. If already bool, just return. If
    not, try to convert and return as a numpy bool array.

    Parameters
    ----------
    column : pandas.Series
        column from dataframe that should be boolean
    column_name : str
        name of column (for error message)

    Returns
    -------
    column : numpy.array
        boolean numpy array.
    """

    # Do a pass trying to infer the datatype of the column. (This is useful if
    # we dropped empty rows that made the original pandas read this column in
    # as a mix of bool and object).
    column = column.infer_objects()

    # If it's not a boolean column, try to turn into one
    if not np.dtype(column.dtypes) is np.dtype(bool):

        # Base message. If everything works great, let user know what
        # happened as warning. If things go awry, use as start of error
        # message
        w = "\n\n"
        w += f"The '{column_name}' column must be boolean (True/False). pandas\n"
        w += "did not recognize the column as boolean, so we're parsing it\n"
        w += "manually by looking for 0/1, yes/no, true/false, etc.\n\n"

        new_column = []
        look_for_true = re.compile("[1yt]",re.IGNORECASE)
        look_for_false = re.compile("[0nf]",re.IGNORECASE)
        for k in column:
            if type(k) is bool:
                is_true = True and k
                is_false = not is_true
                looks_like_a = "bool"
            elif type(k) is str:
                is_true = look_for_true.search(k) is not None
                is_false = look_for_false.search(k) is not None
                looks_like_a = "string"
            elif type(k) is int:
                is_true = (k != 0)
                is_false = (k == 0)
                looks_like_a = "int"
            elif type(k) is float:
                is_true = np.logical_not(np.isclose(k,0))
                is_false = np.isclose(k,0)
                looks_like_a = "float"
            else:
                w += f"Could not figure out how to parse '{k}'\n\n"
                raise ValueError(w)

            if (is_true and is_false) or (not is_true and not is_false):
                w += f"Trying to parse '{k}' as a {looks_like_a}, but\n"
                w += f"could not figure out whether true or false.\n\n"
                raise ValueError(w)
            else:
                new_column.append(is_true)

        # Record newly boolean-ized values
        column = np.array(new_column,dtype=bool)

        # Let user know we manually parsed the keep column...
        print(w)

    return column
