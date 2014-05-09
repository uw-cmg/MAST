#!/bin/bash

# declare a name for this job (replace <jobname> with something more descriptive)
#$ -N ?mast_name?

# request the queue for this job
# replace <queue_name> with morgan.q, morganshort.q, izabela.q
#$ -q ?mast_queue?

# request computational resources for this job as follows
# mvapich2 is current parallel environment. Do not change unless you use a different MPI
# <num> specifies how many processors in total to request. Set it as follows:
#    Replace <num> below with 24 for morgan.q, morganshort.q
#    Replace <num> below with 24 for izabela.q
#$ -pe mvapich2 ?mast_processors?

# request 48 hours of wall time
# for morgan.q, max is 168 hrs; from morganshort.q, max is 24 hrs
#$ -l h_rt=?mast_walltime?

# run the job from the directory of submission.Uncomment only if you don't want the defults.
#$ -cwd
# combine SGE standard output and error files
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
# transfer all your environment variables. Uncomment only if you don't want the defults
#$ -V

# The following is for reporting only. It is not really needed to run the job
# It will show how many processors you get in your output file
echo "Got $NSLOTS processors."


# Use full pathname to make sure we are using the right mpi
#If you are not using the general purpose mpiexec, make sure your mpi environment is properly set up such
#that the correct mpirun is found (you should use the mpirun provided with the compiler
#used to compile the program you are running).

MPI_HOME=/usr/local/mvapich2/intel/1.8.1/bin

$MPI_HOME/mpiexec -n $NSLOTS ?mast_exec?
