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
import numpy as np
import logging
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Lattice
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import dirutil
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
from MAST.ingredients.checker import VaspNEBChecker
from MAST.ingredients.checker import VaspChecker

class WriteIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)

    def no_setup(self):
        """No setup is needed."""
        return

    
    def write_neb(self):
        """Get the parent structures, sort and match atoms, and interpolate.
            Write images to the appropriate folders.
        """
        parentstructures = self.get_parent_structures()
        parentimagestructures = self.get_parent_image_structures()
        if len(parentimagestructures) == 0:
            sxtend = StructureExtensions(struc_work1=parentstructures[0], struc_work2=parentstructures[1])
            image_structures = sxtend.do_interpolation(self.keywords['program_keys']['images'])
        else:
            image_structures = list()
            image_structures.append(parentstructures[0])
            image_structures.extend(parentimagestructures)
            image_structures.append(parentstructures[1])
        if image_structures == None:
            raise MASTError(self.__class__.__name__,"Bad number of images")
        self.checker.keywords['program_keys']['image_structures']=image_structures
        self.checker.set_up_program_input()
        self.place_parent_energy_files()
        self.write_submit_script()
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
        struct_init = self.checker.get_structure_from_file(pfpath_init)
        struct_fin = self.checker.get_structure_from_file(pfpath_fin)
        base_struct = self.keywords['structure']
        mylabel = BaseIngredient.get_my_label(self, "neb_label")
        neblines = self.keywords['program_keys']['neblines'][mylabel]
        sxtendi = StructureExtensions(struc_work1 = struct_init, struc_init = base_struct)
        sorted_init = sxtendi.sort_structure_and_neb_lines(neblines, '00', self.keywords['program_keys']['images']) 
        sxtendf = StructureExtensions(struc_work1 = struct_fin, struc_init = base_struct)
        sorted_fin = sxtendf.sort_structure_and_neb_lines(neblines, str(self.keywords['program_keys']['images'] + 1).zfill(2), self.keywords['program_keys']['images'])
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
                struct_im = self.checker.get_structure_from_file(pfpath)
                base_struct = self.keywords['structure']
                mylabel = BaseIngredient.get_my_label(self, "neb_label")
                neblines = self.keywords['program_keys']['neblines'][mylabel]
                sxtend = StructureExtensions(struc_work1 = struct_im, struc_init = base_struct)
                sorted_im = sxtend.sort_structure_and_neb_lines(neblines, str(imct).zfill(2), self.keywords['program_keys']['images']) 
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
        pffile1.to_file(self.checker.get_path_to_write_neb_parent_energy(1)) #MASTFile contains directory locking.
        pffile2.to_file(self.checker.get_path_to_write_neb_parent_energy(2))
        return



    def write_neb_subfolders(self):
        """Write the static runs to each subfolder.
        """
        myname=self.keywords['name']
        mystr=self.keywords['structure']
        numim = int(self.keywords['program_keys']['images'])
        singlechecker = self.checker
        if self.program == 'vasp':
            nebchecker = VaspNEBChecker(name=self.checker.keywords['name'], program_keys = dict(self.checker.keywords['program_keys']), structure = self.checker.keywords['structure'].copy())
        else:
            raise MASTError(self.__class__.__name__, "NEB checker not supported for program %s") % self.program
        self.checker = nebchecker
        self.write_neb() #Write the POSCAR and OSZICAR files from existing parent-given structures
        self.checker = singlechecker
        for imct in range(1, numim+1):
            subdir = str(imct).zfill(2)
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            self.checker.set_up_program_input()
            self.keywords['name'] = newname
            self.write_submit_script()
        self.keywords['name'] = myname
        self.checker.keywords['name']=myname
        self.checker.keywords['structure']=mystr
        return
    def write_singlerun(self):
        self.checker.set_up_program_input()
        self.write_submit_script()
    def write_singlerun_automesh(self):
        self.checker.scale_mesh(1000)
        self.checker.set_up_program_input()
        self.write_submit_script()
    def write_phonon_multiple(self):
        """Write the multiple phonon files, one for each atom and each direction.
        """
        self.checker.set_up_program_input()
        self.write_submit_script()
        mystructure = self.checker.get_initial_structure_from_directory()
        [pcs,pcr] = self.get_my_phonon_params()
        sxtend = StructureExtensions(struc_work1 = mystructure)
        sdarrlist = sxtend.get_multiple_sd_array(pcs, pcr)
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
            self.checker.keywords['name']=newname
            self.checker.keywords['structure']=mystructure
            self.checker.set_up_program_input()
            self.checker.add_selective_dynamics_to_structure_file(sdarr)
            self.keywords['name'] = newname
            self.write_submit_script()
            self.checker.keywords['name']=myname
            self.checker.softlink_charge_density_file(newname)
            self.checker.softlink_wavefunction_file(newname)
            sct = sct + 1
        self.checker.keywords['name'] = myname
        self.keywords['name']=myname
        

    def write_phonon_single(self):
        """Write the phonon files to a directory.
        """
        self.checker.set_up_program_input()
        self.write_submit_script()
        mystructure = self.checker.get_initial_structure_from_directory()
        [pcs,pcr] = self.get_my_phonon_params()
        sxtend = StructureExtensions(struc_work1 = mystructure)
        sdarr = sxtend.get_sd_array(pcs, pcr)
        if sdarr == None:
            return
        self.checker.add_selective_dynamics_to_structure_file(sdarr)

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
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)
    def ready_singlerun(self):
        return BaseIngredient.is_ready_to_run(self)
    def ready_structure(self):
        if self.directory_is_locked():
            return False
        return self.checker.has_starting_structure_file()

    def ready_defect(self):
        self.logger.warning("ready_defect is deprecated. use ready_structure instead")
        return self.ready_structure()

    def ready_neb_subfolders(self):
        """Make sure all subfolders are ready to run."""
        myname=self.keywords['name']
        notready=0
        numim = int(self.keywords['program_keys']['images'])
        for imct in range(0,numim+2):
            subdir = str(imct).zfill(2)
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            else:
                if dirutil.directory_is_locked(newname):
                    notready = notready + 1
                if not self.checker.is_ready_to_run():
                    notready = notready + 1
            imct = imct + 1
        self.checker.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

    def ready_subfolders(self):
        """Make sure all subfolders are ready to run."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        notready=0
        if len(subdirs) == 0:
            return False
        for subdir in subdirs:
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            if dirutil.directory_is_locked(newname):
                notready = notready + 1
            if not self.checker.is_ready_to_run():
                notready = notready + 1
        self.checker.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

class RunIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)
    def run_singlerun(self, mode='serial'):
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

    def run_defect(self):
        try:
            base_structure = self.checker.get_initial_structure_from_directory() 
        except: #no initial structure
            base_structure = self.keywords['structure'].copy()
            self.logger.warning("No parent structure detected for induce defect ingredient %s. Using initial structure of the recipe." % self.keywords['name'])

        defect_label = BaseIngredient.get_my_label(self, "defect_label")
        defect = self.keywords['program_keys'][defect_label]
        for key in defect:
            if 'subdefect' in key:
                subdefect = defect[key]
                sxtend = StructureExtensions(struc_work1=base_structure)
                base_structure = sxtend.induce_defect(subdefect, defect['coord_type'], defect['threshold'])
            else:
                pass
        self.checker.write_final_structure_file(base_structure)
        return

    def run_strain(self):
        """Strain the lattice.
            Args:
                Looks for mast_strain in input file, for 
                percent strains.
                mast_strain 1.01 1.01 -0.98
            Returns:
                Creates structure file in directory 
        """
        mystructure = self.checker.get_initial_structure_from_directory()
        mystrain = self.keywords['program_keys']['mast_strain']
        sxtend = StructureExtensions(struc_work1 = mystructure)
        strained_structure = sxtend.strain_lattice(mystrain)
        self.checker.write_final_structure_file(strained_structure)
        return


    def run_scale_defect(self):
        try:
            base_structure = self.checker.get_initial_structure_from_directory() 
        except: #no initial structure
            base_structure = self.keywords['structure'].copy()
            self.logger.warning("No parent structure detected for induce defect ingredient %s. Using initial structure of the recipe." % self.keywords['name'])
        scalextend = StructureExtensions(struc_work1=base_structure)
        if not 'mast_scale' in self.keywords['program_keys'].keys():
            raise MASTError(self.__class__.__name__,"No mast_scale ingredient keyword for scaling ingredient %s." % self.keywords['name'])
        scaled = scalextend.scale_structure(int(self.keywords['program_keys']['mast_scale']))
        defect_label = BaseIngredient.get_my_label(self, "defect_label")
        defect = self.keywords['program_keys'][defect_label]
        for key in defect:
            if 'subdefect' in key:
                subdefect = defect[key]
                sxtend = StructureExtensions(struc_work1=scaled, struc_work2=base_structure)
                scaled = sxtend.scale_defect(subdefect, defect['coord_type'], defect['threshold'])
            else:
                pass
        self.checker.write_final_structure_file(scaled)
        return

class IsCompleteIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)
    def complete_structure(self):
        if self.directory_is_locked():
            return False
        return self.checker.has_ending_structure_file()
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
            self.checker.keywords['name']=newname
            self.errhandler.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            elif not self.is_complete():
                notready = notready + 1
            imct = imct + 1
        self.keywords['name']=myname
        self.checker.keywords['name']=myname
        self.errhandler.keywords['name']=myname
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
            self.checker.keywords['name']=newname
            self.errhandler.keywords['name']=newname
            if not self.is_complete():
                notready = notready + 1
        self.keywords['name']=myname
        self.checker.keywords['name']=myname
        self.errhandler.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

class UpdateChildrenIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)
    def _fullpath_childname(self, childname):
        """Make sure the childname has a full path.
            Args:
                childname <str>: child name, like defect_1
            Returns:
                fullpathchild <str>: full path to child, e.g.
                $MAST_SCRATCH/recipefolder/defect_1
        """
        fullpathchild=""
        mydirname = os.path.dirname(self.keywords['name'])
        fullpathchild = os.path.join(mydirname, childname)
        return fullpathchild

    def give_structure(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)

    def give_neb_structures_to_neb(self, childname):
        """Update to ANOTHER NEB."""
        childname = self._fullpath_childname(childname)
        myct=1
        while myct <= self.keywords['program_keys']['images']:
            imno = str(myct).zfill(2)
            impath = os.path.join(self.keywords['name'], imno)
            self.checker.keywords['name'] = impath
            self.checker.forward_final_structure_file(childname,"parent_structure_" + BaseIngredient.get_my_label(self, "neb_label") + '_' + imno)
            myct = myct + 1
    






    def give_saddle_structure(self, childname):
        """Forward the middle image structure."""
        childname = self._fullpath_childname(childname)
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        highenergyimct=0
        imct=0
        highenergy=-1000000000000.0
        #numim = int(self.keywords['program_keys']['images'])
        #middleim = int(numim/2)+1 #returns 1 for 1, 2 for 3 im, 3 for 5 im, etc
        subdirs.sort()
        if self.program == 'vasp_neb' or self.program == 'vasp':
            singlechecker = VaspChecker(name=self.checker.keywords['name'],program_keys = dict(self.checker.keywords['program_keys']),structure = self.checker.keywords['structure'].copy())
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported" % self.program)
        for subdir in subdirs:
            singlechecker.keywords['name'] = subdir
            myenergy = singlechecker.get_energy_from_energy_file()
            if myenergy > highenergy:
                highenergyimct = imct
                highenergy = myenergy
            imct = imct + 1
        imno = str(highenergyimct).zfill(2)
        impath = os.path.join(self.keywords['name'], imno)
        self.checker.keywords['name'] = impath
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
        self.checker.softlink_wavefunction_file(childname)
        return
    def give_phonon_multiple_forces_and_displacements(self,childname):
        self.checker.combine_dynamical_matrix_files(self.keywords['name'])
        self.checker.combine_displacement_files(self.keywords['name'])
        self.give_phonon_single_forces_and_displacements(childname)
    def give_phonon_single_forces_and_displacements(self, childname):
        #Do NOT forward the CONTCAR structure, since the ending CONTCAR contains a displacement in it. Instead, forward the POSCAR
        childname = self._fullpath_childname(childname)
        self.checker.forward_dynamical_matrix_file(childname)
        self.checker.forward_displacement_file(childname)
        self.checker.forward_initial_structure_file(childname, "POSCAR_prePHON")

    def give_structure_and_energy_to_neb(self, childname):
        childname = self._fullpath_childname(childname)
        label = BaseIngredient.get_my_label(self, "defect_label")
        self.checker.forward_final_structure_file(childname,"parent_structure_" + label)
        self.checker.forward_energy_file(childname, "parent_energy_" + label)
    def give_structure_and_restart_files_softlinks(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
        self.checker.softlink_wavefunction_file(childname)
    
    def give_structure_and_restart_files(self, childname):
        self.give_structure_and_restart_files_softlinks(childname)

    def give_structure_and_restart_files_full_copies(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_charge_density_file(childname)
        self.checker.forward_wavefunction_file(childname)
   
    def give_structure_and_charge_density_full_copy(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_charge_density_file(childname)
    
    def give_structure_and_wavefunction_full_copy(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_wavefunction_file(childname)

    def give_structure_and_charge_density_softlink(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
    
    def give_structure_and_wavefunction_softlink(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_wavefunction_file(childname)

