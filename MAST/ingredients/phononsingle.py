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
        Optimize.update_children()
        BaseIngredient.forward_parent_dynmat(self)

    def get_my_label(self, signal="phonon"):
        """Get the calculation label.
            There should be nothing after the label (e.g. _phonon_label, 
            with nothing after it.)
            Args:
                signal <str>: string to signal the beginning of the label.
            Returns:
                The next "_"-delimited piece after the first occurrence of
                signal.
        """
        bname = os.path.basename(self.keywords['name'])
        namesplit = bname.split('_')
        if not (signal in namesplit):
            raise MASTError(self.__class__.__name__, "Phonon label could not be found after signal '%s_'" % signal)
        sidx = namesplit.index(signal)
        label = '_'.join(namesplit[sidx+1:])
        if label == "":
            raise MASTError(self.__class__.__name__, "No label found after signal '%s_' in ingredient name %s" % (signal, bname)) 
        self.label = label
        return label


    def write_files(self):
        """Write the single phonon files. 
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
        if not self.label in self.keywords['program_keys']['phonon'].keys():
            raise MASTError(self.__class__.__name__, "Label %s for phonons not found in phonon input dict." % self.label)

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
        
        pcscoord = np.array(self.phonon_center_site.strip().split(), float)
        pcsarr = pymatgen.util.coord_utils.find_in_coord_list(mystruc.frac_coords, pcscoord,tol)
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
