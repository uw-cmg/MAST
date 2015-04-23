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

def save_old_status_files():
    """Save off old status files from $MAST_SCRATCH
       into $MAST_CONTROL
    """
    mastcontrol=dirutil.get_mast_control_path()
    if not os.path.exists(os.path.join(mastcontrol,"statusfiles")):
        os.mkdir(os.path.join(mastcontrol,"statusfiles"))
    mastscratch=dirutil.get_mast_scratch_path()
    recipedirs=dirutil.immediate_subdirs(mastscratch)
    for recipedir in recipedirs:
        statusfile = MASTFile(os.path.join(mastscratch,recipedir,"status.txt"))
        trydir = os.path.join(mastcontrol,"statusfiles",recipedir)
        if not os.path.exists(trydir):
            os.mkdir(trydir)
        statusfile.to_file("%s/status.txt" % trydir)
    return True

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
            totarfile.data.append(os.path.join("MAST_2/SCRATCH",recipe,ingred)+"\n")
            totarfile.data.append(os.path.join("MAST_2/SCRATCH",recipe,"status.txt") + "\n")
        elif ":archived" in ritem:
            recipe=ritem.split(":")[0]
            totarfile.data.append(os.path.join("MAST_2/ARCHIVE",recipe)+"\n") 
    totarfile.data.append("MAST_2/CONTROL\n")
    totarfile.to_file(os.getcwd() + "/tarthis")
    return

if __name__ == "__main__":
    main()
