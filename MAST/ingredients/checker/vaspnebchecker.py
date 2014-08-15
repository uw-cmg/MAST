##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.cifio import CifParser
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
import os
import shutil
import logging
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import BaseChecker
class VaspNEBChecker(VaspChecker):
    """VASP checker functions
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)
        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])

    def get_path_to_write_neb_parent_energy(self, parent):
        """Get the path into which to write the NEB parent
            energy.
            For VASP, the paths are the 00 and 0(N+1) 
            directories, where N is the number of images.
            Args:
                parent (if integer) <int>: 1 for initial endpoint; 2 for
                              final endpoint
                parent (if string) <str>: destination folder specified as 
                              a string, e.g. "03"
        """
        myname = self.keywords['name']
        myimages = self.keywords['program_keys']['mast_neb_settings']['images']
        if parent == 1:
            return os.path.join(myname, "00", "OSZICAR")
        elif parent == 2:
            return os.path.join(myname, str(int(myimages)+1).zfill(2),"OSZICAR")
        elif len(parent) > 1:
            return os.path.join(myname, parent, "OSZICAR")
        else:
            raise MASTError(self.__class__.__name__,"Parent not specified correctly.")

    def set_up_neb_folders(self, image_structures):
        """Set up NEB folders.
            Args:
               image_structures <list of Structure>: List
                   of image structures
        """
        imct=0
        myname = self.keywords['name']
        if 'mast_coordinates' in self.keywords['program_keys'].keys():
            coordstrucs=self.get_coordinates_only_structure_from_input()
            newstrucs=list()
            sidx = 0 #ex. coordstrucs 0, 1, 2 for 3 images
            while sidx < self.keywords['program_keys']['mast_neb_settings']['images']:
                sxtend = StructureExtensions(struc_work1=image_structures[sidx+1].copy(), name=self.keywords['name'])
                newstrucs.append(sxtend.graft_coordinates_onto_structure(coordstrucs[sidx]))
                sidx = sidx + 1
        while imct < len(image_structures):
            imposcar = Poscar(image_structures[imct])
            num_str = str(imct).zfill(2)
            impath = os.path.join(myname, num_str)
            impospath = os.path.join(myname, "POSCAR_" + num_str)
            if 'mast_coordinates' in self.keywords['program_keys'].keys():
                if imct == 0: #skip endpoint
                    pass 
                elif imct == len(image_structures)-1: #skip other endpt
                    pass
                else:
                    imposcar.structure=newstrucs[imct-1].copy()
            dirutil.lock_directory(myname)
            imposcar.write_file(impospath)
            dirutil.unlock_directory(myname)
            try:
                os.makedirs(impath)
            except OSError:
                self.logger.warning("Directory at %s already exists." % impath)
                return None
            dirutil.lock_directory(impath)
            imposcar.write_file(os.path.join(impath, "POSCAR"))
            dirutil.unlock_directory(impath)
            imct = imct + 1
        return
        
    def is_frozen(self):
        """Check if single VASP non-NEB calculation is frozen.
        """
        return BaseChecker.is_frozen(self, "vasprun.xml")
    def is_complete(self):
        """Check if all images in a VASP NEB calculation are complete.
        """
        dirname = self.keywords['name']
        numim = int(self.keywords['program_keys']['mast_neb_settings']['images'])
        imct=1
        while imct <= numim:
            num_str = str(imct).zfill(2)
            impath = os.path.join(dirname, num_str)
            singlechecker=VaspChecker(name=impath,program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            if not singlechecker.is_complete():
                return False
            imct = imct + 1
        return True

    def is_ready_to_run(self):
        """Check if single VASP NEB ingredient is 
            ready to run.
        """
        dirname = self.keywords['name']
        notready=0
        if not(os.path.isfile(dirname + "/KPOINTS")):
            notready = notready + 1
        if not(os.path.isfile(dirname + "/POTCAR")):
            notready = notready + 1
        if not(os.path.isfile(dirname + "/INCAR")):
            notready = notready + 1
        else: #INCAR exists. Now check POSCARs.
            subdirs = dirutil.walkdirs(dirname)
            if not (len(subdirs) == self.keywords['program_keys']['mast_neb_settings']['images'] + 2):
                notready = notready + 1 #bad number of subdirs
            subdirs.sort()
            if len(subdirs) == 0: #This is an NEB without folders yet
                notready = notready + 1
            else:
                for subdir in subdirs:
                    if not (os.path.isfile(subdir + "/POSCAR")):
                        notready = notready + 1
                if not os.path.isfile(subdirs[0] + "/OSZICAR"):
                    notready = notready + 1
                if not os.path.isfile(subdirs[-1] + "/OSZICAR"):
                    notready = notready + 1
        if not(os.path.isfile(dirname + "/submit.sh")):
            notready = notready + 1
        if notready > 0:
            return False
        else:
            return True
    
    def _vasp_neb_incar_modify(self):
        """Modify the INCAR to add the IMAGES tag back in.
        """
        name=self.keywords['name']
        images = str(self.keywords['program_keys']['mast_neb_settings']['images'])
        my_incar = MASTFile(self.keywords['name'] + "/INCAR")
        my_incar.data.append("IMAGES=%s\n" % images)
        my_incar.to_file(self.keywords['name'] + "/INCAR")
        return my_incar

    def set_up_program_input(self):
        image_structures = self.keywords['program_keys']['mast_neb_settings']['image_structures']
        self.set_up_neb_folders(image_structures)
        self._vasp_kpoints_setup()
        mypotcar = self._vasp_potcar_setup(Poscar(image_structures[0]))
        self._vasp_incar_setup(mypotcar, Poscar(image_structures[0]))
        self._vasp_neb_incar_modify()
        return

    def get_energy_from_energy_file(self):
        """Get the energy from the energy file.
            For VASP neb, this is the last E0 energy from each
            OSZICAR,
            and this function should be used if vasprun.xml
            is corrupted or not available.
            Args:
                mydir <str>: Directory in which to look.
            Returns:
                <str>: all last E0 energies from OSZICAR files.
        """
        myct=0
        mystr=""
        while myct <= self.keywords['program_keys']['mast_neb_settings']['images']+1:
            fullpath = os.path.join(self.keywords['name'], str(myct).zfill(2), "OSZICAR")
            if not os.path.isfile(fullpath):
                raise MASTError(self.__class__.__name__, "No OSZICAR file at %s" % fullpath)
            myosz = MASTFile(fullpath)
            mye0 = myosz.get_segment_from_last_line_match("E0", "E0=","d E =")
            mye0 = float(mye0)
            mystr = mystr + "%3.3f" % mye0 + ';'
            myct = myct + 1
        mystr=mystr[0:-1] #remove last semicolon
        return mystr

    def is_started(self):
        """See if the ingredient has been started on
            the queue.
        """
        if os.path.isfile(os.path.join(self.keywords['name'],'01','OUTCAR')):
            return True
        else:
            return False
