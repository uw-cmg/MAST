import numpy as np

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
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def write_files(self):
        """Write the single phonon file. 
        """
        Optimize.write_files(self)
        mystructure = self.get_structure_from_directory(self.keywords['name'])
        sdarr = self.get_sd_array(mystructure)
        if sdarr == None:
            return
        self.add_selective_dynamics_to_structure(self, sdarr)

    def get_neighbor_array(self, mystruc, tol=1e-1):
        """
            Get a neighbor-index array.
            Use program_keywords 'phonon_center_site' and 
            'phonon_neighbor_radius' to limit the number of phonons calculated.
            self.program_keys['phonon_center_site'] should be a list of 
                coordinates.
                    If the key is missing, all atoms will be taken into account.
                    If there is more than one atom in the list, all of these
                        atoms which are found will be taken into account.
            self.program_keys['phonon_neighbor_radius'] 
                    should be a positive float.
                    If the key is missing or 0, nothing extra happens.
                    If the key is present and nonzero, then all atoms in a
                        radius around EACH site found in phonon_center_site
                        will also be taken into account.
            Args:
                mystruc <Structure>: pymatgen Structure
                tol <float>: Tolerance for match-searching.
        """
        if not 'phonon_center_site' in self.program_keys.keys():
            return None
        pcstotarr=None
        pcsarr=None
        for pcs in self.program_keys['phonon_center_site']:
            pcsarr = pymatgen.util.coord_utils.find_in_coord_list(mystruc.frac_coords, np.array(pcs.strip().split()),tol)
            if pcstotarr == None:
                pcstotarr = pcsarr
            else:
                pcstotarr = np.concatenate([pcstotarr, pcsarr])
        
        uniqsites = np.unique(pcstotarr)
        
        if len(uniqsites) == 0:
            raise MASTError(self.__class__.__name__, "No sites found for phonon centering.")

        if not 'phonon_neighbor_radius' in self.program_keys.keys():
            return uniqsites
        
        nrad = self.program_keys['phonon_neighbor_radius']
        if nrad == 0:
            return uniqsites
        if nrad < 0:
            raise MASTError(self.__class__.__name__, "Phonon neighbor radius should not be less than zero!")

        nbtotarr=None
        for pcs in uniqsites:
            neighbors = mystruc.get_neighbors(mystruc[pcs], nrad, True)
            if nbtotarr == None:
                nbtotarr = neighbors
            else:
                np.concatenate([nbtotarr, neighbors])
        alltotarr = np.concatenate([uniqsites, nbtotarr])
        allsites = np.unique(alltotarr)
        return alltotarr

    def get_sd_array(self, mystruc):
        """Create a selective dynamics array.
            Args:
                mystruc <Structure>: pymatgen Structure
        """
        if not 'phonon_center_site' in self.program_keys.keys():
            return None
        mynbarr = self.get_neighbor_array(mystruc)
        mysd = np.zeros([mystruc.num_sites,3],bool)
        for myn in mynbarr:
            mysd[myn]=np.ones(3,bool)
        return mysd
