__description__ = \
"""
Shared functions for running pipelines.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

import topiary

def _run_and_print(function,kwargs,step_counter,out_file_string,human_string):
    """
    Wrapper function that calls another function in a stereotyped way.

    Parameters
    ----------
        function: function to call. code assumes this function returns a
                  topiary dataframe.
        kwargs: keyword arguments to pass to the function
        step_counter: number output file with this counter
        out_file_string: write the resulting dataframe to a csv file with
                        the format step_counter_out_file_string-datframe.csv
        human_string: Write this to stdout for human readability.

    Return
    ------
        dataframe output of function, incremented step counter
    """

    print("-------------------------------------------------------------------")
    print(human_string)
    print("-------------------------------------------------------------------")
    print("",flush=True)

    df = function(**kwargs)
    topiary.write_dataframe(df,f"{step_counter:02d}_{out_file_string}-dataframe.csv")

    print("",flush=True)

    return df, step_counter + 1
