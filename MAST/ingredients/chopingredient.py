############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os
import shutil
import numpy as np

import pymatgen
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import dirutil
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.pmgextend import structure_extensions
class WriteIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)


    def write_defect(self):
        work_dir = '/'.join(self.keywords['name'].split('/')[:-1])
        name = self.keywords['name'].split('/')[-1]
        print "write_files:", name
        BaseIngredient.get_structure_from_directory(self, self.keywords['name']) #TTM+1 changed to use BaseIngredient method

        defect_label = BaseIngredient.get_my_label(self, "defect_label")
        #print 'GRJ DEBUG: defect_label =', defect_label
        #defect_label = 'defect_' + name.split('/')[-1].split('_')[-1]
        #print 'GRJ DEBUG: defect_label (regex) =', defect_label

        #self.metafile.write_data('debug', [name, name.split('/')])
        defect = self.keywords['program_keys'][defect_label]
        #print 'Defect in write_files:', defect

        base_structure = self.keywords['structure'].copy()
        for key in defect:
            if 'subdefect' in key:
                subdefect = defect[key]
                base_structure = structure_extensions.induce_defect(base_structure, subdefect, defect['coord_type'], defect['threshold'])
                #base_structure = modified_structure
            else:
                pass

        if self.keywords['program'].lower() == 'vasp':
            #myposcar = Poscar(modified_structure)
            myposcar = Poscar(base_structure)
            #print "poscar OK"
            self.lock_directory()
            #print "lock OK"
            myposcar.write_file('%s/%s/CONTCAR' % (work_dir, name))
            #print "Write sucessful"
            self.unlock_directory()
            #print "Unlock sucessful"
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])

        return
    
    def write_neb(self):
        """Get the parent structures, sort and match atoms, and interpolate.
            Write images to the appropriate folders.
        """
        parentstructures = self.get_parent_structures()
        parentimagestructures = self.get_parent_image_structures()
        if len(parentimagestructures) == 0:
            image_structures = structure_extensions.do_interpolation(parentstructures, self.keywords['program_keys']['images'])
        else:
            image_structures = list()
            image_structures.append(parentstructures[0])
            image_structures.extend(parentimagestructures)
            image_structures.append(parentstructures[1])
        if image_structures == None:
            raise MASTError(self.__class__.__name__,"Bad number of images")
        if self.keywords['program'] == 'vasp':
            BaseIngredient.set_up_program_input_neb(self, image_structures)
            self.place_parent_energy_files()
            self.write_submit_script()
        else:
            raise MASTError(self.__class__.__name__,"Program not supported. No setup accomplished.")
        return





        
    def get_parent_structures(self):
        """Assume that parents have written files
            named 'parent_structure_<N>'. 
            For VASP these are CONTCAR-type files.
            Returns:
                [struct_init, struct_fin]: pymatgen Structure objects
        """
        header = os.path.join(self.keywords['name'], "parent_structure_")
        mylabel = BaseIngredient.get_my_label(self, "neb_label").split("-")
        pfpath_init = header + mylabel[0]
        pfpath_fin = header + mylabel[1]
        if not os.path.isfile(pfpath_init):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_init)
        if not os.path.isfile(pfpath_fin):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_fin)
        struct_init = BaseIngredient.get_structure_from_file(self, pfpath_init)
        struct_fin = BaseIngredient.get_structure_from_file(self, pfpath_fin)
        base_struct = self.keywords['structure']
        mylabel = BaseIngredient.get_my_label(self, "neb_label")
        neblines = self.keywords['program_keys']['neblines'][mylabel]
        sorted_init = structure_extensions.sort_structure_and_neb_lines(struct_init, base_struct, neblines, 0) 
        sorted_fin = structure_extensions.sort_structure_and_neb_lines(struct_fin, base_struct, neblines, 1)
        return [sorted_init, sorted_fin]

    def get_parent_image_structures(self):
        """A low-mesh NEB may have written files
            named 'parent_structure_<N-N>_0N'. 
            For VASP these are CONTCAR-type files.
            Returns:
                list of <Structure>: list of pymatgen Structure objects
        """
        header = "parent_structure_"
        numim = self.keywords['program_keys']['images']
        imct = 1
        imstrs=list()
        while imct <= numim:
            pfpath=""
            for myfile in os.listdir(self.keywords['name']):
                if (header in myfile) and (str(imct).zfill(2) in myfile):
                    pfpath = os.path.join(self.keywords['name'],myfile)
            if pfpath == "":
                pass
            else:
                struct_im = BaseIngredient.get_structure_from_file(self, pfpath)
                base_struct = self.keywords['structure']
                mylabel = BaseIngredient.get_my_label(self, "neb_label")
                neblines = self.keywords['program_keys']['neblines'][mylabel]
                sorted_im = structure_extensions.sort_structure_and_neb_lines(struct_im, base_struct, neblines, 0) 
                imstrs.append(sorted_im)
            imct = imct + 1
        if len(imstrs) > 0 and not (len(imstrs) == numim):
            raise MASTError(self.__class__.__name__, "Incomplete number of forwared images found!")
        return imstrs

    def place_parent_energy_files(self):
        """Assume parents have written files parent_energy_<N>.
            Copy these files into the 00 and 0N directories.
        """
        header = os.path.join(self.keywords['name'], "parent_energy_")
        mylabel=BaseIngredient.get_my_label(self, "neb_label").split("-")
        pfpath1= header + mylabel[0]
        pfpath2= header + mylabel[1]
        pffile1=MASTFile(pfpath1)
        pffile2=MASTFile(pfpath2)
        pffile1.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self,1)) #MASTFile contains directory locking.
        pffile2.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, 2))
        return



    def write_neb_subfolders(self):
        """Write the static runs to each subfolder.
        """
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
        self.write_neb() #Write the POSCAR files from existing parent-given structures
        return
    def write_singlerun(self):
        self.set_up_program_input()
        self.write_submit_script()
    def write_phonon_multiple(self):
        """Write the multiple phonon files, one for each atom and each direction.
        """
        self.set_up_program_input()
        self.write_submit_script()
        mystructure = BaseIngredient.get_structure_from_directory(self, self.keywords['name'])
        [pcs,pcr] = self.get_my_phonon_params()
        sdarrlist = structure_extensions.get_multiple_sd_array(pcs, pcr, mystructure)
        if sdarrlist == None:
            raise MASTError(self.__class__.__name__, "No phonons to run!")
        sct=1
        myname=self.keywords['name']
        for sdarr in sdarrlist:
            newname = os.path.join(myname,"phon_%s" % str(sct).zfill(2))
            try:
                os.mkdir(newname)
            except OSError:
                pass
            self.keywords['name']=newname
            self.keywords['structure']=mystructure
            self.set_up_program_input()
            self.add_selective_dynamics_to_structure(sdarr)
            self.write_submit_script()
            self.forward_extra_restart_files(myname, newname)
            sct = sct + 1
        self.keywords['name']=myname
        

    def write_phonon_single(self):
        """Write the phonon files to a directory.
        """
        self.set_up_program_input()
        self.write_submit_script()
        mystructure = BaseIngredient.get_structure_from_directory(self, self.keywords['name'])
        [pcs,pcr] = self.get_my_phonon_params()
        sdarr = structure_extensions.get_sd_array(pcs, pcr, mystructure)
        if sdarr == None:
            return
        self.add_selective_dynamics_to_structure(sdarr)

    def get_my_phonon_params(self):
        """Get phonon parameters from 
            ['program_keys']['phonon'][label]['phonon_center_site'] and
            ['program_keys']['phonon'][label]['phonon_center_radius'] 
            Returns:
                [phonon_center_site, phonon_center_radius]
        """
        if not 'phonon' in self.keywords['program_keys'].keys():
            return [None, None]
        mylabel = BaseIngredient.get_my_label(self, "phonon_label")
        if not mylabel in self.keywords['program_keys']['phonon'].keys():
            raise MASTError(self.__class__.__name__, "Label %s for phonons not found in phonon input dict for %s" % (mylabel, self.keywords['name']))

        myphdict = dict(self.keywords['program_keys']['phonon'][mylabel])
        if not 'phonon_center_site' in myphdict.keys():
            return [None, None]
        phonon_center_site = myphdict['phonon_center_site']
        if not 'phonon_center_radius' in myphdict.keys():
            phonon_center_radius = None
        else:
            phonon_center_radius = myphdict['phonon_center_radius']
        return [phonon_center_site,phonon_center_radius]


class IsReadyToRunIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
    def ready_singlerun(self):
        return BaseIngredient.is_ready_to_run(self)
    def ready_defect(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'].lower() == 'vasp':
            if os.path.exists(self.keywords['name'] +'/POSCAR'):
                return True
            else:
                return False
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])

    def ready_neb_subfolders(self):
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

    def ready_subfolders(self):
        """Make sure all phonon subfolders are ready to run."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        notready=0
        if len(phondirs) == 0:
            return False
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            if not BaseIngredient.is_ready_to_run(self):
                notready = notready + 1
        self.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

class RunIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
    def run_noqueue(self, mode='noqsub'):
        self.write_files()
        return True
    def run_singlerun(self, mode='serial', curdir=os.getcwd()):
        return BaseIngredient.run(self, mode)
    def run_neb_subfolders(self):
        """Run all image subfolders."""
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
    def run_subfolders(self):
        """Run all subfolders."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            BaseIngredient.run(self,"serial")
        self.keywords['name']=myname
        return



class IsCompleteIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
    def complete_structure(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'].lower() == 'vasp':
            if os.path.exists(self.keywords['name'] +'/CONTCAR'):
                return True
            else:
                return False
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])
    def complete_singlerun(self):
        return BaseIngredient.is_complete(self)
    def complete_neb_subfolders(self):
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

    def complete_subfolders(self):
        """Make sure all subfolders are complete."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        notready=0
        if len(phondirs) == 0:
            return False
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            if not BaseIngredient.is_complete(self):
                notready = notready + 1
        self.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

class UpdateChildrenIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
    def give_structure(self, childname):
        self.forward_parent_structure(self.keywords['name'], childname)

    def give_neb_structures_to_neb(self, childname):
        """Update to ANOTHER NEB."""
        myct=1
        while myct <= self.keywords['program_keys']['images']:
            imno = str(myct).zfill(2)
            impath = os.path.join(self.keywords['name'], imno)
            self.forward_parent_structure(impath, childname,"parent_structure_" + BaseIngredient.get_my_label(self, "neb_label") + '_' + imno)
            myct = myct + 1
    






    def give_saddle_structure(self, childname):
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
                self.forward_parent_structure(newname, childname)
                self.forward_extra_restart_files(newname, childname)
            imct = imct + 1
        return
    def give_phonon_multiple_forces_and_displacements(self,childname):
        self.combine_dynmats()
        shutil.copy(os.path.join(self.keywords['name'],"DYNMAT_combined"),
            os.path.join(self.keywords['name'],"DYNMAT"))
        self.combine_displacements()
        shutil.copy(os.path.join(self.keywords['name'],"XDATCAR_combined"),
            os.path.join(self.keywords['name'],"XDATCAR"))
        self.give_phonon_single_forces_and_displacements(childname)
    def give_phonon_single_forces_and_displacements(self, childname):
        #Do NOT forward the CONTCAR structure, since the ending CONTCAR contains a displacement in it. Instead, forward the POSCAR
        self.forward_parent_dynmat(self.keywords['name'], childname)
        self.forward_parent_initial_structure(self.keywords['name'],childname, "POSCAR_prePHON")

    def give_structure_and_energy_to_neb(self, childname):
        label = BaseIngredient.get_my_label(self, "defect_label")
        self.forward_parent_structure(self.keywords['name'], childname,"parent_structure_" + label)
        self.forward_parent_energy(self.keywords['name'], childname, "parent_energy_" + label)
    def give_structure_and_restart_files(self, childname):
        self.forward_parent_structure(self.keywords['name'], childname)
        self.forward_extra_restart_files(self.keywords['name'], childname)
   

