#!/bin/env python
############
# TTM 2014-08-28 nshd_mast
############
import os
import sys
import optparse
from MAST.utility.mastfile import MASTFile
from MAST.utility import dirutil
from MAST.recipe.recipeplan import RecipePlan

def main():
    """Compares status files"""
    parser=optparse.OptionParser()
    parser.add_option('-m','--mode',dest='mode',default='None',
                       help='save or compare status files')
    (myopt, myarg) = parser.parse_args()
    if myopt.mode == "save":
        save_old_status_files()
    elif myopt.mode == "compare":
        compare_old_and_new_status_files()
    elif myopt.mode == "zip":
        zip_only_changed()
    else:
        print "No option selected. Exiting."
    return


def save_recipe_files(listfilename="",whichfilename=""):
    """Save off old files from recipe in $MAST_SCRATCH
       into $MAST_CONTROL
       Args:
        listfilename <str>: File name for the list of files to be saved
        whichfilename <str>: Short file name to be saved (identical file
                                for all recipes)
    """
    if listfilename == "":
        raise NameError("No list file name given.")
    if whichfilename == "":
        raise NameError("No file name given for which files to save.")
    mastcontrol=dirutil.get_mast_control_path()
    if not os.path.exists(os.path.join(mastcontrol,listfilename)):
        os.mkdir(os.path.join(mastcontrol,listfilename))
    mastscratch=dirutil.get_mast_scratch_path()
    recipedirs=dirutil.immediate_subdirs(mastscratch)
    for recipedir in recipedirs:
        statusfile = MASTFile(os.path.join(mastscratch,recipedir,whichfilename))
        trydir = os.path.join(mastcontrol,listfilename,recipedir)
        if not os.path.exists(trydir):
            os.mkdir(trydir)
        statusfile.to_file("%s/%s" % (trydir, whichfilename))
    return True

def save_ingredient_files(listfilename="",whichfilename=""):
    """Save off old files from ingredients in $MAST_SCRATCH
       into $MAST_CONTROL
       Args:
        listfilename <str>: File name for the list of files to be saved
        whichfilename <str>: Short file name to be saved (identical file
                                for all recipes)
    """
    if listfilename == "":
        raise NameError("No list file name given.")
    if whichfilename == "":
        raise NameError("No file name given for which files to save.")
    mastcontrol=dirutil.get_mast_control_path()
    if not os.path.exists(os.path.join(mastcontrol,listfilename)):
        os.mkdir(os.path.join(mastcontrol,listfilename))
    mastscratch=dirutil.get_mast_scratch_path()
    recipedirs=dirutil.immediate_subdirs(mastscratch)
    for recipedir in recipedirs:
        recipefulldir=os.path.join(mastscratch,recipedir)
        ingreddirs=dirutil.immediate_subdirs(recipefulldir)
        for ingreddir in ingreddirs:
            statfilepath = os.path.join(recipefulldir, ingreddir, whichfilename)
            if os.path.exists(statfilepath):
                statusfile = MASTFile(statfilepath)
                trydir = os.path.join(mastcontrol,listfilename,recipedir,ingreddir)
                if not os.path.exists(trydir):
                    os.mkdir(trydir)
                statusfile.to_file("%s/%s" % (trydir, whichfilename))
    return True

def save_old_status_files():
    """Save off old status files from $MAST_SCRATCH
       into $MAST_CONTROL
    """
    return save_recipe_files("statusfiles","status.txt")

def save_old_change_status_files():
    """Save off old change_status.txt files from ingredient directories
        into $MAST_CONTROL
    """
    return save_ingredient_files("changestatusfiles","change_status.txt")

def compare_old_and_new_status_files():
    """Compare old and new status files from $MAST_CONTROL
        and $MAST_SCRATCH.
        Write an indication whether the status files have changed,
        not changed. Use "archived" if a recipe has been archived or is
        missing from SCRATCH for some other reason.
    """
    rdict=dict()
    mastcontrol=dirutil.get_mast_control_path()
    mastscratch=dirutil.get_mast_scratch_path()
    recipedirs=dirutil.immediate_subdirs(os.path.join(mastcontrol,"statusfiles"))
    for recipedir in recipedirs:
        mystatus="unknown"
        rdict[recipedir]=dict()
        changelist=list()
        if not os.path.exists(os.path.join(mastscratch,recipedir)):
            mystatus="archived"
	else:
            scratchstatusfile = MASTFile(os.path.join(mastscratch,recipedir,"status.txt"))
            controlstatusfile = MASTFile(os.path.join(mastcontrol,"statusfiles",recipedir,"status.txt"))
            if scratchstatusfile.data == controlstatusfile.data:
                mystatus="unchanged"
            else:
                mystatus="changed"
                myidx=0
                while myidx < len(scratchstatusfile.data):
                    oldline = controlstatusfile.data[myidx]
                    newline = scratchstatusfile.data[myidx]
                    if "#" in oldline:
                        pass
                    else:
                        ingred = oldline.split(":")[0].strip()
                        oldstatus = oldline.split(":")[1].strip()
                        newstatus = newline.split(":")[1].strip()
                        if (oldstatus == "P") and (newstatus == "P"):
                             rdict[recipedir][ingred]="AVOID"
                        elif (oldstatus == "C") and (newstatus == "C"):
                             rdict[recipedir][ingred]="AVOID"
                        else:
                             rdict[recipedir][ingred]="send"
                    myidx = myidx + 1
        rdict[recipedir]["MAIN"]=mystatus
    recipedirs=dirutil.immedate_subdirs(os.path.join(mastcontrol,"changestatusfiles"))
    for recipedir in recipedirs:
        ingreddirs=dirutil.immediate_subdirs(os.path.join(mastcontrol,"changestatusfiles",recipedir))
        for ingreddir in ingreddirs:
            oldchangestatus=MASTFile(os.path.join(mastcontrol,"changestatusfiles",recipedir,ingreddir,"change_status.txt"))
            newchangestatus=MASTFile(os.path.join(mastscratch,recipedir,ingreddir,"change_status.txt"))
            if not (oldchangestatus.data == newchangestatus.data):
                rdict[recipedir]["MAIN"]="changed"
                rdict[recipedir][ingreddir]="send"
    rchangefile = MASTFile()
    for rdir in rdict.keys():
        for key, value in rdict[rdir].iteritems():
            rchangefile.data.append("%s:%s:%s\n" % (rdir,key,value))
    rchangefile.to_file(os.path.join(mastcontrol,"recipechangefile.txt"))
    return True

def zip_only_changed():
    """Zip only the changed recipes"""
    mastcontrol=dirutil.get_mast_control_path()
    mastscratch=dirutil.get_mast_scratch_path()
    rchangefile = MASTFile(os.path.join(mastcontrol,"recipechangefile.txt"))
    #import shutil
    #shutil.rmtree(os.path.join(mastcontrol,"statusfiles"))
    totarfile = MASTFile()
    for ritem in rchangefile.data:
        if ":send" in ritem:
            recipe=ritem.split(":")[0]
            ingred=ritem.split(":")[1]
            totarfile.data.append(os.path.join("MAST/SCRATCH",recipe,ingred)+"\n")
            totarfile.data.append(os.path.join("MAST/SCRATCH",recipe,"status.txt") + "\n")
        elif ":archived" in ritem:
            recipe=ritem.split(":")[0]
            totarfile.data.append(os.path.join("MAST/ARCHIVE",recipe)+"\n") 
    totarfile.data.append("MAST/CONTROL\n")
    totarfile.to_file(os.getcwd() + "/tarthis")
    return

if __name__ == "__main__":
    main()
