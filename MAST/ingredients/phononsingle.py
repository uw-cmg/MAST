import numpy as np
import pymatgen
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.utility.mastobj import MASTObj
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.pmgextend import vasp_extensions
from MAST.ingredients.optimize import Optimize
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata

import os
import shutil
#import pdb
#TTM

class PhononSingle(Optimize):
    """
        Attributes:
            self.phonon_center_site <np.array of float>: coords of center site
            self.phonon_center_radius <float>: radius around center site
            self.label <str>: label of this calculation
    """
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
    def update_children(self):
        #Do NOT forward the structure, since the ending CONTCAR contains a displacement in it.
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_dynmat(self.keywords['name'], childname)

    def get_my_label(self):
        """Get the phonon label from the metadata file.
            Returns:
                plabel <str>: phonon label
        """
        myname = self.keywords['name']
        mymeta = Metadata(metafile=os.path.join(myname, "metadata.txt"))
        mylabel = mymeta.search_data("phononlabel")
        if mylabel == "":
            raise MASTError(self.__class__.__name__,
                "No metadata file for tag phononlabel.")
        plabel = mylabel[1]
        self.label=plabel
        return 
        

    def write_files(self):
        """Write the multiple phonon files to different directories
        """
        self.get_my_label()
        self.get_my_phonon_params()
        self.set_up_program_input()
        self.write_submit_script()
        mystructure = self.get_structure_from_directory(self.keywords['name'])
        sdarr = self.get_sd_array(mystructure)
        if sdarr == None:
            return
        self.add_selective_dynamics_to_structure(sdarr)

    def get_my_phonon_params(self):
        """Get phonon parameters from 
            ['program_keys']['phonon'][label]['phonon_center_site'] and
            ['program_keys']['phonon'][label]['phonon_center_radius'] 
            and set them into class attributes.
        """
        if not 'phonon' in self.keywords['program_keys'].keys():
            self.phonon_center_site = None
            self.phonon_center_radius = None
            return
        self.get_my_label()
        if not self.label in self.keywords['program_keys']['phonon'].keys():
            raise MASTError(self.__class__.__name__, "Label %s for phonons not found in phonon input dict for %s" % (self.label, self.keywords['name']))

        myphdict = dict(self.keywords['program_keys']['phonon'][self.label])
        if not 'phonon_center_site' in myphdict.keys():
            self.phonon_center_site = None
            self.phonon_center_radius = None
            return
        self.phonon_center_site = myphdict['phonon_center_site']
        if not 'phonon_center_radius' in myphdict.keys():
            self.phonon_center_radius = None
        else:
            self.phonon_center_radius = myphdict['phonon_center_radius']
        return

    def get_neighbor_array(self, mystruc, tol=1e-1):
        """
            Get a neighbor-index array.
            Use program_keywords 'phonon_center_site' and 
            'phonon_center_radius' to limit the number of phonons calculated.
            ['program_keys']['phonon'][label]['phonon_center_site'] 
                    should be a coordinate
                    If the key is missing, all atoms will be taken into account.
            ['program_keys']['phonon'][label]['phonon_center_radius'] 
                    should be a positive float in ANGSTROMS (Not fractional.)
                    If the key is missing or 0, nothing extra happens.
                    If the key is present and nonzero, then all atoms in a
                        radius around EACH site found in phonon_center_site
                        will also be taken into account.
            Args:
                mystruc <Structure>: pymatgen Structure
                tol <float>: Tolerance for match-searching.
        """
        if self.phonon_center_site == None:
            return None
        print "TTM DEBUG: centersite: ", self.phonon_center_site.strip().split()
        print "TTM DEBUG: MYSTRUC: ", mystruc 
        pcscoord = np.array(self.phonon_center_site.strip().split(), float)
        pcsarr = pymatgen.util.coord_utils.find_in_coord_list(mystruc.frac_coords, pcscoord,tol)
        print "TTM DEBUG PCSarr: ", pcsarr
        uniqsites = np.unique(pcsarr)
        
        if len(uniqsites) == 0:
            raise MASTError(self.__class__.__name__, "No sites found for phonon centering.")

        if self.phonon_center_radius == None:
            return uniqsites
        
        nrad = float(self.phonon_center_radius)
        if nrad == 0:
            return uniqsites
        if nrad < 0:
            raise MASTError(self.__class__.__name__, "Phonon center radius should not be less than zero!")

        nbtotarr=None
        for pcs in uniqsites:
            neighbors = mystruc.get_neighbors(mystruc[pcs], nrad, True)
            if nbtotarr == None:
                nbtotarr = neighbors
            else:
                np.concatenate([nbtotarr, neighbors])
        nbsitelist=list()
        for nbr in nbtotarr:
            nbsitelist.append(nbr[-1])
        nbsitelist = np.array(nbsitelist)
        alltotarr = np.concatenate([uniqsites, nbsitelist])
        allsites = np.unique(alltotarr)
        return allsites

    def get_sd_array(self, mystruc):
        """Create a selective dynamics array.
            Args:
                mystruc <Structure>: pymatgen Structure
        """
        if self.phonon_center_site == None:
            return None
        mynbarr = self.get_neighbor_array(mystruc)
        mysd = np.zeros([mystruc.num_sites,3],bool)
        for myn in mynbarr:
            mysd[myn]=np.ones(3,bool)
        return mysd
