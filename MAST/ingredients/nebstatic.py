import numpy as np
import pymatgen
from MAST.ingredients.neb import NEB
from MAST.ingredients.baseingredient import BaseIngredient
import os
import shutil
from MAST.utility import dirutil

class NEBStatic(NEB):
    """
        Perform static calculation on each image.
    """
    def __init__(self, **kwargs):
        NEB.__init__(self, **kwargs)

    def update_children(self):
        """Forward the middle image structure."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        imct = 0
        numim = int(self.keywords['program_keys']['images'])
        middleim = int(numim/2)+1 #returns 1 for 1, 2 for 3 im, 3 for 5 im, etc
        for subdir in subdirs:
            if not (imct == middleim):
                pass
            else:
                newname = os.path.join(myname, subdir)
                for childname in self.keywords['child_dict'].iterkeys():
                    self.forward_parent_structure(newname, childname)
            imct = imct + 1
        return

    def write_files(self):
        """Write the static runs to each subfolder.
        """
        NEB.write_files(self) #Write the POSCAR files
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        imct = 0
        numim = int(self.keywords['program_keys']['images'])
        for subdir in subdirs:
            if imct == 0 or imct > numim:
                pass
            else:
                newname = os.path.join(myname, subdir)
                try:
                    os.mkdir(newname)
                except OSError:
                    pass
                self.keywords['name']=newname
                self.set_up_program_input()
                self.write_submit_script()
            imct = imct + 1
        self.keywords['name']=myname
        return
    def is_complete(self):
        """Make sure all subfolders are complete."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        if len(subdirs) == 0:
            return False
        notready=0
        imct = 0
        numim = int(self.keywords['program_keys']['images'])
        for subdir in subdirs:
            newname = os.path.join(myname, subdir)
            self.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            elif not BaseIngredient.is_complete(self):
                #print "TTM DEBUG: Not ready at imct %1i" % imct
                notready = notready + 1
            imct = imct + 1
        self.keywords['name']=myname
        #print "TTM DEBUG: Notready: ", notready
        if notready == 0:
            return True
        else:
            return False

    def is_ready_to_run(self):
        """Make sure all subfolders are ready to run."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        if len(subdirs) == 0:
            return False
        notready=0
        imct = 0
        numim = int(self.keywords['program_keys']['images'])
        for subdir in subdirs:
            newname = os.path.join(myname, subdir)
            self.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            elif not BaseIngredient.is_ready_to_run(self):
                notready = notready + 1
            imct = imct + 1
        self.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

    def run(self):
        """Run all subfolders."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        imct = 0
        numim = int(self.keywords['program_keys']['images'])
        for subdir in subdirs:
            if imct == 0 or imct > numim:
                pass
            else:
                newname = os.path.join(myname, subdir)
                self.keywords['name']=newname
                BaseIngredient.run(self,"serial")
            imct = imct + 1
        self.keywords['name']=myname
        return

