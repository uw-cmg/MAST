#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
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

def direct_shell_command():
    return "//share/apps/vasp5.2_cNEB"

def queue_submission_command():
    return "qsub"

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
            'Q': Queued
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
     
    qmat = queuetext.split()
    #These are approximate translations of the actual
    #    job status. 
    #The point is to prevent MAST from resubmitting a job
    #    which is already on the queue. The user should check
    #    qstat manually in order to see the actual job status.
    qpc = qmat[4].lower()
    if 'e' in qpc: #Some error found
        return "E" 
    elif 'd' in qpc: #In the process of being deleted
        return "R"
    elif 'h' in qpc: #Held somehow
        return "H"
    elif 't' in qpc: #Transferring
        return "H"
    elif 's' in qpc: #Suspended
        return "H"
    elif 'q' in qpc: #Queued
        return "Q"
    elif 'r' in qpc: #Running
        return "R"
    else: #Not sure what this string is
        return "E"

def extract_submitted_jobid(string):
    """
        Extract the job ID from a string returned when the job is submitted
        INPUTS:
            string <str> = string with job ID somewhere in it
        OUTPUTS:
            <int> = job ID as integer
    """
    #Looks like: Your job 46745 ("jobname") has been submitted
    if string == "":
        return None
    return int(string.split()[2])

def queue_snap_command():
    """
        Return the command for producing a queue snapshot.
        INPUTS:
            None
        OUTPUTS:
            Command for producing a queue snapshot
    """
    uname = os.getenv("USER")
    return "qstat -u %s" % uname

def get_approx_job_error_file(jobid):
    """
        Return the approximate job error file name
            Args:
                jobid <int>: Job ID
    """
    if jobid == None:
        return ".e"
    else:
        return ".e%s" % str(jobid)
