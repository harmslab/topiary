#!/usr/bin/env python
"""
Run sphinx-apidocs and clean up. This should be run locally and the resulting
rst files committed to github for conversion to proper docs on readthedocs.
"""

import os
import glob
import re

def main():
    """
    Main function. Build .rst from source code and clean up.
    """

    # Run sphinx-apidoc on
    os.system("sphinx-apidoc -f ../../topiary -o .")

    # This wackiness cleans up the sphinx-apidoc so it's much cleaner and easier
    # to read
    topiary_rst = glob.glob("topiary.*.rst")

    for rst in topiary_rst:

        header = []
        body = []
        footer = []

        wipe_out_next = False
        get_footer = False
        counter = 0
        with open(rst) as f:
            for line in f:

                if counter < 3:
                    header.append(line)
                    counter += 1
                    continue

                if line.startswith("Submodules") or line.startswith("Subpackages"):
                    wipe_out_next = True
                    continue

                if wipe_out_next:
                    wipe_out_next = False
                    continue

                if line.startswith("Module contents"):
                    get_footer = True

                if get_footer:
                    footer.append(line)
                    continue

                line = re.sub(" module","",line)
                body.append(line)


        header[0] = re.sub(" package","",header[0])

        final_out = header[:]
        final_out.append("\n")
        final_out.extend(footer[2:])
        final_out.append("\n")
        final_out.extend(body)
        final_out.append("\n")

        g = open(rst,"w")
        g.writelines(final_out)
        g.close()

if __name__ == "__main__":
    main()
