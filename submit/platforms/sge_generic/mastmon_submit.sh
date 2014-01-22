#!/bin/bash

# declare a name for this job (replace <jobname> with something more descriptive)
#$ -N mastmon

# request the queue for this job
#$ -q default.q

# request 48 hours of wall time
#$ -l h_rt=04:00:00

# Specifies the interpreting shell for this job to be the Bash shell
#$ -S /bin/bash

# run the job from the directory of submission.Uncomment only if you don't want the defults.
#$ -cwd
# combine SGE standard output and error files
#$ -o mastmon_submission_output
#$ -e mastmon_submission_errors
# transfer all your environment variables. Uncomment only if you don't want the defults
#$ -V

# Put the command to run one per line
python runmast.py >> $MAST_CONTROL/mastoutput 2> $MAST_CONTROL/errormast
