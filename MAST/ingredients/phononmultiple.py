import numpy as np
import pymatgen
from MAST.ingredients.phononsingle import PhononSingle
import os
class PhononMultiple(PhononSingle):
    """
        Split out phonons into individual atoms and directions.
        Attributes:
            self.phonon_center_site <np.array of float>: coords of center site
            self.phonon_center_radius <float>: radius around center site
            self.label <str>: label of this calculation
    """
    def __init__(self, **kwargs):
        PhononSingle.__init__(self, **kwargs)
    def update_children(self):
        #Do NOT forward the structure, since the ending CONTCAR contains a displacement in it.
        self.recombine_dynmat()
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_dynmat(self.keywords['name'], childname)

    def write_files(self):
        """Write the multiple phonon files, one for each atom and each direction.
        """
        self.get_my_label()
        self.get_my_phonon_params()
        self.set_up_program_input()
        self.write_submit_script()
        mystructure = self.get_structure_from_directory(self.keywords['name'])
        sdarrlist = self.get_multiple_sd_array(mystructure)
        if sdarrlist == None:
            raise MASTError(self.__class__.__name__, "No phonons to run!")
        sct=1
        myname=self.keywords['name']
        for sdarr in sdarrlist:
            newname = os.path.join(myname,"phon_%s" % str(sct).zfill(2))
            os.mkdir(newname)
            self.keywords['name']=newname
            self.keywords['structure']=mystructure
            self.set_up_program_input()
            self.add_selective_dynamics_to_structure(sdarr)
            self.write_submit_script()
            sct = sct + 1
        self.keywords['name']=myname

    def is_ready_to_run(self):
        """Make sure all subfolders are ready to run."""
        print "Make sure subfolders are ready to run."
        return False

    def run(self):
        """Run all subfolders."""
        print "Run all subfolders."
        return False

    def get_multiple_sd_array(self, mystruc):
        """Create a selective dynamics array.
            Args:
                mystruc <Structure>: pymatgen Structure
            Returns:
                mysdlist <list>: list of SD arrays
        """
        if self.phonon_center_site == None:
            return None
        mynbarr = self.get_neighbor_array(mystruc)
        mysdlist=list()
        for myn in mynbarr:
            for myct in range(0,3):
                mysd = np.zeros([mystruc.num_sites,3],bool)
                mysd[myn]=np.zeros(3,bool)
                mysd[myn][myct]=1
                mysdlist.append(mysd)
        return mysdlist
