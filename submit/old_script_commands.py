#!/usr/bin/env python
# Purpose: Modify and get information out of submission scripts.
# Authors: Tam Mayeshiba
# 11/9/12 created
import os
import sys
import math #TA to calc num_nodes
import MAST.utility.fileutil

system = "bardeen" #dlx, bardeen, curie, nersc, ranger

def get_num_nodes(scriptpath):
    """
        Get the number of nodes from the script
        INPUTS:
            scriptpath <str> = full path to the script file
        OUTPUTS:
            numnodes <int> = number of nodes
    """
    script=fileutil.MASTFile(scriptpath)
    if script.data == []:
        print "Error getting number of nodes: Script empty at ", scriptpath
        return None

    if system == "bardeen" or system == "curie":###################TA 11/14/12: BARDEEN/CURIE CASE##############
        node_num=script.get_line_match_number("nodes=")
        node_line = script.data[node_num-1]

        index = node_line.find("nodes=") #TA 11/14/12 find the line number with the num_nodes in it
        index2 = node_line.find(":") #TA find where ends in case of double digits
        numnodes_str = node_line[index+6:index2] #TA add length of search string (6 for bardeen)
    
        return int(numnodes_str)
    
    if system == "nersc": ###################TA 11/14/12: NERSC CASE################
        num_nodes_set = [4,6,9,8,12,12,12,12,12] #TA deciding to just store the set of nodes used-too complic to rig. calc
        num_proc_sets = get_num_procs(scriptpath)/24 #Theoretical number of  nodes - want more like sqrt num procs though

        if num_proc_sets > 9:
            return 12
    
        return num_nodes_set[num_proc_sets-1]

    if system == "dlx":###################TTM 20121206 UKY DLX##################### 
        node_num=script.get_line_match_number("#SBATCH -N")
        node_line = script.data[node_num-1]
        return int(node_line.split()[2])

def get_num_procs(scriptpath):
    """
        Get the number of processors form the script
        INPUTS:
            scriptpath <str> = full path to the script file
        OUTPUTS:
            numprocs <int> = number of processors
    """
    script=fileutil.MASTFile(scriptpath)
    if script.data == []:
        print "Error getting number of processors: Script empty at ", scriptpath
        return None

    if system == "bardeen" or system == "curie":###################TA 11/14/12: BARDEEN CASE#############
        proc_line_num=script.get_line_match_number("ppn=")
        proc_line = script.data[proc_line_num-1]

        index = proc_line.find("ppn=") #TA 11/14/12 find the line number with the num_procs in it
        index2 = proc_line.find(",") #TA find where ends
        numprocs_str = proc_line[index+4:index2] #TA add length of search string
    
        return int(numprocs_str)

    if system == "nersc": ###################TA: NERSC CASE##############
        proc_line_num=script.get_line_match_number("mppwidth=")
        proc_line = script.data[proc_line_num-1]

        index = proc_line.find("mppwidth=") #TA 11/14/12 find the line number with the num procs in it

        numprocs_str = proc_line[index+9:] #TA add length of search string (9 to end for nersc)

        return int(numprocs_str)

    if system == "dlx": #############TTM: UKY DLX ##################
        proc_line_num=script.get_line_match_number("#SBATCH -n")
        # This is actually the number of procs per node (16) * num nodes
        proc_line = script.data[proc_line_num-1]
        return int(proc_line.split()[2])

def get_memory(scriptpath):
    """
        Get the number of processors form the script
        INPUTS:
            scriptpath <str> = full path to the script file
        OUTPUTS:
            memory <str> = memory string, where applicable
    """
    print "TESTING"
    return "2000mb"

def modify_jobname(scriptpath, jobname):
    """
        Change the jobname in the submission script.
        INPUTS:
            scriptpath <str> = full path to submission script file
            jobname <str> = new job name to use
        OUTPUTS:
            Modifies the file at scriptpath to use the new jobname.
    """
    if not '.sh' in scriptpath:
        print "Error: ",scriptpath, "does not point to a .sh file."
        return None
    myscript=fileutil.MASTFile(scriptpath)
    if myscript.data == []:
        print "Error modifying script at ", scriptpath
        return None

    if system == "bardeen": #### BARDEEN / CURIE / NERSC #########
        namenum=myscript.get_line_match_number("PBS -N")
        myscript.modify_file_by_line_number(namenum, "R", "#PBS -N " + jobname + "\n")

    if system == "dlx": #### UKY DLX ###########
        namenum=myscript.get_line_match_number("#SBATCH -J")
        myscript.modify_file_by_line_number(namenum, "R", "#SBATCH -J " + jobname + "\n")
    
    myscript.to_file(scriptpath)
    return True

def multiply_nodes_and_procs(scriptpath, multiplier):
    """
        Multiplies the number of nodes by the multiplier.
        Adjusts the number of processors if necessary.
        INPUTS:
            scriptpath <str> = full path to submission script file
            multiplier <str> = digit indicating multiplier to use
        OUTPUTS:
            Modifies the file at scriptpath to use a new number of nodes and processors.
    """
    if not '.sh' in scriptpath:
        print "Error: ",scriptpath, "does not point to a .sh file."
        return None
    myscript=fileutil.MASTFile(scriptpath)
    if myscript.data == []:
        print "Error modifying script at ", scriptpath
        return None
    numnodes = get_num_nodes(scriptpath)*int(multiplier)
    numprocs = get_num_procs(scriptpath)
    memstr = get_memory(scriptpath)
    
    if system == "bardeen":##### BARDEEN ONLY ######
        if numnodes > 2:
            numprocs = 12 #switch to morgan2

            namenum=myscript.get_line_match_number("PBS -l nodes")
            myscript.modify_file_by_line_number(namenum, "R", "#PBS -l nodes=" + str(numnodes) + ':ppn=' + str(numprocs) + ',pvmem='+  memstr + "\n")
            myscript.to_file(scriptpath)
    
        if numnodes > 2:
            modify_queue(scriptpath, 'morgan2')
            modify_walltime(scriptpath, '96')
        return True

    if system == "dlx":##### UKY DLX ######
        numprocs = numprocs*int(multiplier)
        namenum=myscript.get_line_match_number("#SBATCH -N")
        myscript.modify_file_by_line_number(namenum, "R", "#SBATCH -N " + str(numnodes) + "\n")
        namenum=myscript.get_line_match_number("#SBATCH -n")
        myscript.modify_file_by_line_number(namenum, "R", "#SBATCH -n " + str(numprocs) + "\n")
        myscript.to_file(scriptpath)
        return True

    if system == "nersc": ###################TA: NERSC CASE### 
        ############just multiply num procs by multiplier
        if not '.sh' in scriptpath:
            print "Error: ",scriptpath, "does not point to a .sh file."
            return None
        myscript=fileutil.MASTFile(scriptpath)
        if myscript.data == []:
            print "Error modifying script at ", scriptpath
            return None
        numprocs = get_num_procs(scriptpath)*multiplier

        namenum=myscript.get_line_match_number("PBS -l mppwidth")
        myscript.modify_file_by_line_number(namenum, "R", "#PBS -l mppwidth=" + str(numprocs) + "\n")
        myscript.to_file(scriptpath)
    
        return True

def modify_walltime(scriptpath, walltime):
    """
        Modifies the walltime.
        INPUTS:
            scriptpath <str> = full path to submission script file
            walltime <str> = walltime in hours
        OUTPUTS:
            Modifies the file at scriptpath to use a new walltime.
    """
    if not '.sh' in scriptpath:
        print "Error: ",scriptpath, "does not point to a .sh file."
        return None
    myscript=fileutil.MASTFile(scriptpath)
    if myscript.data == []:
        print "Error modifying script at ", scriptpath
        return None

    if system == "bardeen" or system == "curie" or system == "nersc": ####### BARDEEN /CURIE /NERSC ############
        linenum=myscript.get_line_match_number("PBS -l walltime")
        myscript.modify_file_by_line_number(linenum, "R", "#PBS -l walltime=" + str(walltime) + ':00:00' + "\n")
        myscript.to_file(scriptpath)
        return True

    if system == "dlx":########UKY DLX ##############
        linenum=myscript.get_line_match_number("#SBATCH -t")
        myscript.modify_file_by_line_number(linenum, "R", "#SBATCH -t " + str(walltime) + ':00:00' + "\n")
        myscript.to_file(scriptpath)
        return True

def modify_queue(scriptpath, queue):
    """
        Modifies the queue
        INPUTS:
            scriptpath <str> = full path to submission script file
            queue <str> = name of queue to use
        OUTPUTS:
            Modifies the file at scriptpath to use a new queue.
    """
    if not '.sh' in scriptpath:
        print "Error: ",scriptpath, "does not point to a .sh file."
        return None
    myscript=fileutil.MASTFile(scriptpath)
    if myscript.data == []:
        print "Error modifying script at ", scriptpath
        return None

    if system == "bardeen" or system == "curie" or system == "nersc": ###### BARDEEN/CURIE/NERSC ############
        linenum=myscript.get_line_match_number("PBS -q")
        myscript.modify_file_by_line_number(linenum, "R", "#PBS -q " + str(queue) + "\n")
        myscript.to_file(scriptpath)
        return True

    if system == "dlx":############ UKY DLX ################
        print "ERROR: Check with UKY to see if changing queue is an option."
        return None
        #####Currently not an option? Use queue "Compute"
