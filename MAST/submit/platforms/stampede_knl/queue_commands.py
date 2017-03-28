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
    return "vasp"

def queue_submission_command():
    return "sbatch"

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
    findstr = string.strip()
    return int(findstr.split()[-1])

def queue_snap_command():
    """
        Return the command for producing a queue snapshot.
        INPUTS:
            None
        OUTPUTS:
            Command for producing a queue snapshot
    """
    uname=os.getenv("USER")
    return "squeue -u %s" % uname

def get_approx_job_error_file(jobid):
    """
        Return the approximate job error file name
            Args:
                jobid <int>: Job ID
    """
    if jobid == None:
        return "slurm"
    else:
        return "slurm-%s.out" % str(jobid)
