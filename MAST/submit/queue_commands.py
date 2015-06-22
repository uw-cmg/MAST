#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
# Purpose: Submit a job to the queue. This function is highly platform-specific and may need to be modified.
# 10/25/12 TTM created
# 12/6/12 TTM added UKY DLX commands
import os
import sys
import math #TA to calc num_nodes
import time
import importlib
import subprocess
import shutil
import logging
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import loggerutils

mast_control = dirutil.get_mast_control_path()

myqc = "MAST.submit.platforms.%s.queue_commands" % dirutil.get_mast_platform()
my_queue_commands = importlib.import_module(myqc)

def direct_shell_command(scriptname="submit.sh"):
    myds = my_queue_commands.direct_shell_command()
    return "%s %s" % (myds, scriptname)

def queue_submission_command(scriptname="submit.sh"):
    myqs = my_queue_commands.queue_submission_command()
    return "%s %s" % (myqs, scriptname)

def write_to_submit_list(mydir):
    """Write an entry to the submission list in 
        $MAST_CONTROL/submitlist
        Args:
            mydir <str>: Directory which includes submission
                        script submit.sh for a single 
                        calculation to be submitted to the
                        queue
    """
    mast_control = dirutil.get_mast_control_path() #setting here instead of globally allows tests to run in isolated test_control folder
    submitlist=os.path.join(mast_control, "submitlist")
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
        Adds a job number to the top of the "jobids" file in each
        ingredient directory.
    """
    submitlist=os.path.join(mast_control, "submitlist")
    if not os.path.isfile(submitlist):
        print "No submission list at %s" % submitlist
        return
    submitfile=MASTFile(submitlist)
    subentries=list(submitfile.data)
    subentries.sort()
    submitted=dict()
    subcommand = queue_submission_command()
    for subentry in subentries: # each entry is a directory name
        if len(subentry) == 0:
            continue
        subentry = subentry.strip()
        if len(subentry) == 0:
            continue
        if not os.path.isdir(subentry):
            continue
        elif subentry in submitted.keys():
            continue
        lastjobid = get_last_jobid(subentry)
        jstatus="not yet determined"
        if not (lastjobid == None):
            jstatus = get_job_status_from_queue_snapshot(subentry, lastjobid)
            if jstatus.lower() in ['r','q','h','e']: #running, queued, held, error
                submitted[subentry]="Already on queue with status %s" % jstatus
                continue
        os.chdir(subentry)
        subme=subprocess.Popen(subcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subme.wait()
        status=subme.communicate()[0]
        write_to_jobids_file(subentry, status)
        submitted[subentry]=status
    print_submitted_dict(submitted)
    os.chdir(mast_control)

def write_to_jobids_file(subentry, status):
    """Write the job id to a jobids file in the ingredient directory.
        Args:
            subentry <str>: Folder from which the job was submitted 
                            (ingredient folder)
            status <str>: Job submission response from queueing system
    """
    jobid = extract_submitted_jobid(status)
    jpath = "%s/jobids" % subentry
    if not os.path.isfile(jpath):
        jfile = MASTFile()
    else:
        jfile = MASTFile(jpath)
    jfile.data.insert(0,str(jobid) + '\n')
    jfile.to_file(jpath)
    return

def get_job_status_from_queue_snapshot(ingpath, jobid):
    """Match a jobid to the queue_snapshot
        Args:
            ingpath <str>: ingredient path
            jobid <int>: job id
        Returns:
            Queue status for the job, according to $MAST_CONTROL/queue_snapshot
    """
    myqs = MASTFile("%s/queue_snapshot" % mast_control)
    jobstatus = ""
    for qsline in myqs.data:
        if str(jobid) in qsline:
            jobstatus = queue_status_from_text(jobid, qsline)
    return jobstatus

def get_last_jobid(ingpath):
    """Get the last job id from the jobids file in the ingredient path
        Args:
            ingpath <str>: ingredient path
        Returns:
            jobid <int>: job id, if available
            None if no jobids file was found
    """
    jpath = "%s/jobids" % ingpath
    if not os.path.isfile(jpath):
        return None
    jfile = MASTFile(jpath)
    if len(jfile.data) == 0:
        return None
    if 'None' in jfile.data[0]:
        return None
    return int(jfile.data[0].strip())

def clear_submission_list():
    """Clear all entries from the submission list at
        $MAST_CONTROL/submitlist
    """
    submitlist=os.path.join(mast_control, "submitlist")
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
    subprint=os.path.join(mast_control, "submitted")
    subrecent=os.path.join(mast_control, "just_submitted")
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
    return my_queue_commands.queue_status_from_text(jobid, queuetext)

def extract_submitted_jobid(string):
    """
        Extract the job ID from a string returned when the job is submitted
        INPUTS:
            string <str> = string with job ID somewhere in it
        OUTPUTS:
            <int> = job ID as integer
    """
    return my_queue_commands.extract_submitted_jobid(string)

def queue_snap_command():
    """
        Return the command for producing a queue snapshot.
        INPUTS:
            None
        OUTPUTS:
            Command for producing a queue snapshot
    """
    return my_queue_commands.queue_snap_command()

def get_job_error_file(ingpath):
    """
        Return the job error file full path
            Args:
                ingpath <str>: ingredient path
    """
    logger = loggerutils.get_mast_logger("mast queue commands")
    jobid = get_last_jobid(ingpath)
    tryfile = my_queue_commands.get_approx_job_error_file(jobid)
    #logger.info("Try this search string: %s" % tryfile)
    if not os.path.isdir(ingpath):
        return None
    dircontents = os.listdir(ingpath)
    tryitems = list()
    for diritem in dircontents:
        if tryfile in diritem:
            tryitems.append(os.path.join(ingpath, diritem))
    if len(tryitems) == 0:
        logger.warning("No job error file found for %s" % ingpath)
        return None
    if len(tryitems) > 1:
        logger.warning("More than one job error file found for %s. Using first file %s" % (ingpath, tryitems[0]))
    return tryitems[0]
    

