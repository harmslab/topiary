"""
Script to be executed from the command line that returns all names of processors
allocated after an mpirun call. Useful to get a list of hosts allocated by
SLURM (or other job manager) that can then be sent to parallel mpirun calls
with fewer numbers of allocated slots.
"""

import mpi4py
from mpi4py import MPI

comm = MPI.COMM_WORLD
name = MPI.Get_processor_name()
print(name)
