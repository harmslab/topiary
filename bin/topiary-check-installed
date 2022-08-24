#!/usr/bin/env python3
"""
Check for installed packages on command line.
"""

import topiary

import os

def main(argv=None):

    required = topiary._private.software_requirements

    print(70*"=")
    print("Checking for installed external software")
    print(70*"=")
    print("",flush=True)

    to_check = []
    for r in required:
        to_check.append({"program":r,
                         "min_version":required[r],
                         "must_pass":False})

    try:
        topiary._private.installed.validate_stack(to_check)
    except RuntimeError:
        print("\nNot all external programs are visible to topiary. Please make")
        print("sure they are installed and that the $PATH environment variable")
        print("has the directories containing the installed software.\n")
        print("The current $PATH visible to topiary contains the following")
        print("directories:")

        path_dirs = os.environ["PATH"].split(os.pathsep)
        for d in path_dirs:
            print(f"    {d}")
        print("")

if __name__ == "__main__":
    main()
