__decription__ = \
"""
Light wrapper with muscle, compatible with muscle 5 and 5.1.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-02"

import subprocess

def run_muscle(input_fasta,
               output_fasta,
               muscle_binary="muscle"):
    """
    Run muscle for sequence alignment.

    input_fasta: input fasta file to align
    output_fasta: output fasta file to store alignment
    muscle_binary: location of muscle binary (default assumes 'muscle' command
                   is in the PATH).
    """

    # Build old and new format command lists
    old_cmd = [muscle_binary]
    old_cmd.extend(["-in",input_fasta])
    old_cmd.extend(["-out",output_fasta])

    new_cmd = [muscle_binary]
    new_cmd.extend(["-align",input_fasta])
    new_cmd.extend(["-output",output_fasta])

    success = False
    cmd_to_try = [new_cmd,old_cmd]
    bad_stdout = []
    bad_stderr = []
    for cmd in cmd_to_try:

        ret = subprocess.run(cmd,capture_output=True)

        # Check for error on return
        if ret.returncode != 0:
            bad_stdout.append("".join([line for line in ret.stdout.decode()]))
            bad_stderr.append("".join([line for line in ret.stderr.decode()]))
        else:
            success = True

    if not success:
        err = "\nTried to run muscle in two ways, using the old (3.8) format and\n"
        err += "the new (5.1) format. Both failed.\n\n."
        err += "stdout for new format:\n{bad_stdout[0]}\n\n"
        err += "stderr for new format:\n{bad_stderr[0]}\n\n"
        err += "stdout for old format:\n{bad_stdout[1]}\n\n"
        err += "stderr for old format:\n{bad_stderr[1]}\n\n"
        raise RuntimeError(err)

    print(f"Success. Output written to {output_fasta}")
