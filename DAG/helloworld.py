#!/usr/bin/env python
"""
Parallel Hello World
"""

from mpi4py import MPI
import sys
import time
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

def main(argv):
    if len(argv) <2:
        dur = 5
    dur = int(argv[1]) # duration of each process
    sys.stdout.write(
        "Hello, World! I am process %d of %d on %s.\n"
        % (rank, size, name))
    
    sys.stdout.write(
        "Process %d of %d on %s runs for %d.\n"
        % (rank, size, name, dur))

    time.sleep(dur)
    sys.stdout.write(
        "Process %d of %d on %s ended.\n"
        % (rank, size, name))

if __name__ == "__main__":
    sys.exit(main(sys.argv))
