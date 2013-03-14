#!/usr/bin/python
import sys
import os

def main(argv):
    """Script to submit mpi jobs with different names, programs, and arguments.
        Example :
        ./mpisub.py JOBNAME program argument
        ./mpisub.py JOBNAME python program argument 
    """
    
    jobname = argv[1]
    arg = ""
    ispython = False
    
    if argv[2] !='python':
        prog = argv[2]
        for s in argv[3:]:
            arg += s+" "
    else:
        prog = argv[3]
        for s in argv[4:]:
            arg += s+" "
            

    f = open('mytest.sh','w')
    script = "#!/bin/bash\n\n# declare a name for this job to be sample_job\n"
    script+="#PBS -N %s\n" %jobname

    script+="""# request 4 processes
#PBS -l nodes=4
# request 4 hours of wall time
#PBS -l walltime=04:00:00
# combine PBS standard output and error files
#PBS -j oe
#PBS -k eo
    
#How many procs do I have?
NN=`cat $PBS_NODEFILE | wc -l`
echo \"Processors received = \"$NN
echo \"script running on host `hostname`\"
    
#cd into the directory where I typed qsub
cd $PBS_O_WORKDIR
echo \"PBS_NODEFILE\"
cat $PBS_NODEFILE
    
#Type in commands to run. Replace a.out with the program name
#to run.
/opt/mpiexec/bin/mpiexec -comm pmi """
    if argv[2] == "python":
        script +="  /share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python "
    script += prog+" "+arg+ "\n"
    f.write(script)
    f.close()
    os.system("chmod o+x mytest.sh")
    os.system("qsub mytest.sh")
    os.system("qstat -u hwkim")
if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
