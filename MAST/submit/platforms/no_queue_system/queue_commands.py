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
    return "bash"

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
    if len(qmat) < 5:
        return 'E'
    return qmat[4] #TTM 1/17/12 qstat returns differently than qstat -a does

    #if system == "ranger":####Ranger
    #    ranger_out = queuetext.split()[4]
    #    if (ranger_out == 'wq') or (ranger_out == 'qw'):
    #        final_out = 'Q'
    #    elif ranger_out == 'r':
    #        final_out = 'R'
    #    return final_out

    #if system == "dlx":####UKy DLX
    #    queuetext = queuetext.strip()
    #    dlx_out = queuetext.split()[4] #indexing starts at 0
    #    if (dlx_out == 'PD'):
    #        final_out = 'Q'
    #    elif dlx_out in ['CG','CD','R']:
    #        final_out = 'R'
    #    return final_out

def extract_submitted_jobid(string):
    """
        Extract the job ID from a string returned when the job is submitted
        INPUTS:
            string <str> = string with job ID somewhere in it
        OUTPUTS:
            <int> = job ID as integer
    """
    return None

def queue_snap_command():
    """
        Return the command for producing a queue snapshot.
        INPUTS:
            None
        OUTPUTS:
            Command for producing a queue snapshot
    """
    return ""

def get_approx_job_error_file(jobid):
    """
        Return the approximate job error file name
            Args:
                jobid <int>: Job ID
    """
    return "ingredient_error"
