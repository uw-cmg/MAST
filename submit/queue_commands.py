#!/usr/bin/env python
# Purpose: Submit a job to the queue. This function is highly platform-specific and may need to be modified.
# Authors: Tam Mayeshiba
# 10/25/12 TTM created
# 12/6/12 TTM added UKY DLX commands
import os
import sys
from MAST.utility import dirutil
import subprocess
import time
from MAST.utility import MASTFile
system = "bardeen" #dlx, curie, bardeen, nersc, ranger
compute_node = True

def direct_shell_command():
    return "//share/apps/vasp5.2_cNEB"

def queue_submission_command():
    return "qsub submit.sh"

def write_to_submit_list(mydir):
    """Write an entry to the submission list in 
        $MAST_CONTROL/submitlist
        Args:
            mydir <str>: Directory which includes submission
                        script submit.sh for a single 
                        calculation to be submitted to the
                        queue
    """
    control=dirutil.get_mast_control_path()
    submitlist=os.path.join(control, "submitlist")
    if os.path.isfile(submitlist):
        submitfile=MASTFile(submitlist)
    else:
        submitfile=MASTFile()
    submitfile.data.append(mydir + "\n")
    submitfile.to_file(submitlist)
    return

def submit_from_submission_list():
    """Submit all entries from the submission list at
        $MAST_CONTROL/submitlist
    """
    control=dirutil.get_mast_control_path()
    submitlist=os.path.join(control, "submitlist")
    if not os.path.isfile(submitlist):
        print "No submission list at %s" % submitlist
        return
    submitfile=MASTFile(submitlist)
    subentries=list(submitfile.data)
    subentries.sort()
    submitted=dict()
    subcommand = queue_submission_command()
    for subentry in subentries: # each entry is a directory name
        subentry = subentry.strip()
        if not os.path.isdir(subentry):
            pass
        elif subentry in submitted.keys():
            pass
        else:
            os.chdir(subentry)
            subme=subprocess.Popen(subcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subme.wait()
            status=subme.communicate()[0]
            submitted[subentry]=status
    print_submitted_dict(submitted)
    os.chdir(control)
        
def clear_submission_list():
    """Clear all entries from the submission list at
        $MAST_CONTROL/submitlist
    """
    control=dirutil.get_mast_control_path()
    submitlist=os.path.join(control, "submitlist")
    if not os.path.isfile(submitlist):
        print "No submission list at %s" % submitlist
        return
    submitfile=MASTFile(submitlist)
    submitfile.data=list()
    submitfile.data.append("\n")
    submitfile.to_file(submitlist)

def print_submitted_dict(submitted):
    """Print a dictionary of all runs submitted into
        $MAST_CONTROL/submitted
        Args:
            submitted <dict>: Dictionary of submitted runs,
                with key as the directory name.
    """
    control=dirutil.get_mast_control_path()
    subprint=os.path.join(control, "submitted")
    subrecent=os.path.join(control, "just_submitted")
    if os.path.isfile(subprint):
        subprintfile=MASTFile(subprint)
    else:
        subprintfile=MASTFile()
    subrecentfile=MASTFile()
    keylist=submitted.keys()
    keylist.sort()
    subprintfile.data.append(time.asctime()+ "\n")
    subrecentfile.data.append(time.asctime()+ "\n")
    for key in keylist:
        subprintfile.data.append("%s:%s\n" % (key, submitted[key]))
        subrecentfile.data.append("%s:%s\n" % (key, submitted[key]))
    subprintfile.to_file(subprint)
    subrecentfile.to_file(subrecent)
    return
    


def queue_status_from_text(jobid, queuetext):
    """
        Returns the queue status from a line with the job ID information
        INPUTS:
            jobid <int> = job ID
            qstr <str> = queue snapshot line containing job ID
        OUTPUTS:
            <str> = queue status 
            'E': Error
            'X': Not found
            'R': Running
            'W': Waiting
    """
    if (queuetext == None): 
        return 'X'
    if len(queuetext) == 0:
        return 'X'
    begin = queuetext.find(str(jobid))
    if begin == -1:
        return 'X'
    queuetext = queuetext[begin:]
    end = queuetext.find('\n')
    if end > -1:
        queuetext = queuetext[:end]

    if system == "ranger":####Ranger
        ranger_out = queuetext.split()[4]
        if (ranger_out == 'wq') or (ranger_out == 'qw'):
            final_out = 'Q'
        elif ranger_out == 'r':
            final_out = 'R'
        return final_out

    if system == "bardeen" or system == "curie":###Bardeen, Curie, etc.
        qmat = queuetext.split()
        if len(qmat) < 5:
            return 'E'
        return qmat[4] #TTM 1/17/12 qstat returns differently than qstat -a does
    
    if system == "dlx":####UKy DLX
        queuetext = queuetext.strip()
        dlx_out = queuetext.split()[4] #indexing starts at 0
        if (dlx_out == 'PD'):
            final_out = 'Q'
        elif dlx_out in ['CG','CD','R']:
            final_out = 'R'
        return final_out

def extract_submitted_jobid(string):
    """
        Extract the job ID from a string returned when the job is submitted
        INPUTS:
            string <str> = string with job ID somewhere in it
        OUTPUTS:
            <int> = job ID as integer
    """
    ####Bardeen and Curie:
    if system == "bardeen" or system == "curie":
        final = ""
        for char in string:
            if(char != '.'):
                final += char
            else:
                return int(final)

    if system == "ranger":####Ranger:
        findme=0
        findstr=""
        findme = string.find("Your job") # ..... Your job 12345678 has been submitted.
        if findme == -1: #TTM+1 2/11/12 if not found, return negative 1
            return -1
        findstr = string[findme+9:]
        splitstr=[]
        splitstr = findstr.split()
        return int(splitstr[0])

    if system == "dlx":###UKY DLX:
        findstr = string.strip()
        return int(findstr.split()[3])

def queue_snap_command():
    """
        Return the command for producing a queue snapshot.
        INPUTS:
            None
        OUTPUTS:
            Command for producing a queue snapshot
    """
    qcmd=""
    if(compute_node):
        if(system == "ranger"):####Ranger, from compute node
            #
            qcmd = "ssh -f ranger 'qstat'"

        if(system == "bardeen"):####Bardeen, from compute node
            qcmd = "ssh -f bardeen 'qstat'" 
    
        if(system == "curie"):####Curie, from compute node
            qcmd = "ssh -f curie 'qstat'"
    else:###Direct from headnode
        qcmd = 'qstat 2> /dev/null' #TTM 052212 sometimes requires /bin/bash
    ###Direct from headnode, UKY DLX:
        if(system == "dlx"):
            uname=os.getenv("USER")
            qcmd = "squeue -u " + uname
    return qcmd

