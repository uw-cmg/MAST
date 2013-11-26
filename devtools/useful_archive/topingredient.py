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

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.ingredients.baseingredient import BaseIngredient


class InduceDefect(BaseIngredient):
    def __init__(self, **kwargs):
        #TTM move allowed_keys into __init__ and match order of other sections
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'child_dict': (dict, dict(), 'Dictionary of children'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

        #print 'Initializing InduceDefect'
        #if (self.keywords['coordtype'] == 'cartesian'):
        #    self.keywords['position'] = self._cart2frac(self.keywords['position'])

    def induce_defect(self, base_structure, defect, coord_type, threshold):
        """Creates a defect, and returns the modified structure
            mast_defect is a dictionary like this: 
            'defect_1': {'symbol': 'cr', 'type': 'interstitial', 
                        'coordinates': array([ 0. ,  0.5,  0. ])}}
            'defect_2': {'symbol': 'o', 'type': 'vacancy', 
                        'coordinates': array([ 0.25,  0.25,  0.25])}
            'coord_type': 'fractional' 
        """
        #base_structure = self.get_new_structure()
        #base_structure = self.keywords['structure']
        #print 'Defect in induce_defect', defect
        #print 'base_structure in induce_defect', base_structure
        struct_ed = StructureEditor(base_structure) #should be updated using get_new_structure)
        symbol = defect['symbol'].title() #Cap first letter

        # If we have cartesian coordinates, then we convert them to fractional here.
        if ('cartesian' in coord_type):
            defect['coordinates'] = self._cart2frac(defect['coordinates'])

        if (defect['type'] == 'vacancy'):
            print 'Creating a %s vacancy at %s' % (symbol, str(defect['coordinates']))

            #print defect['coordinates']
            index = find_in_coord_list(base_structure.frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)
            #print base_structure.frac_coords
            #print 'Index of deleted atom is', index
            struct_ed.delete_site(index)
        elif (defect['type'] == 'interstitial'):
            print 'Creating a %s interstitial at %s' % (symbol, str(defect['coordinates']))

            struct_ed.append_site(symbol,
                                  defect['coordinates'],
                                  coords_are_cartesian=False,
                                  validate_proximity=True)
        elif (defect['type'] in ['antisite', 'substitution']):
            print 'Creating a %s antisite at %s' % (symbol, str(defect['coordinates']))

            index = find_in_coord_list(base_structure.frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)

            struct_ed.replace_site(index, symbol)
        else:
            raise RuntimeError('Defect type %s not supported' % defect['type'])

        return struct_ed.modified_structure

    def _cart2frac(self, position):
        """Converts between cartesian coordinates and fractional coordinates"""
        fractional =  self.keywords['structure'].lattice.get_fractional_coords(position)
        for i in range(len(fractional)):
            if (fractional[i] < 0.0):
                fractional[i] += 1.0
            elif (fractional[i] > 1.0):
                fractional[i] -= 1.0
        return fractional

    def write_files(self):
        work_dir = '/'.join(self.keywords['name'].split('/')[:-1])
        name = self.keywords['name'].split('/')[-1]
        print "write_files:", name
        self.get_new_structure()

        defect_label = self.metafile.read_data('defect_label')
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
                base_structure = self.induce_defect(base_structure, subdefect, defect['coord_type'], defect['threshold'])
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
    
    def is_ready_to_run(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'].lower() == 'vasp':
            if os.path.exists(self.keywords['name'] +'/POSCAR'):
                return True
            else:
                return False
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])

    def run(self, mode='noqsub'):
        if self.is_ready_to_run():
            self.write_files()
        return True

    def is_complete(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'].lower() == 'vasp':
            if os.path.exists(self.keywords['name'] +'/CONTCAR'):
                return True
            else:
                return False
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])

    def update_children(self):
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname)

    def get_new_structure(self):
        if self.keywords['program'].lower() == 'vasp':
            if os.path.isfile(self.keywords['name'] + "/POSCAR"):
                myposcar = Poscar.from_file(self.keywords['name'] + "/POSCAR")
                self.keywords['structure'] = myposcar.structure
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported." % self.keywords['program'])
        return
    
    def get_my_number(self):
        """For defect in the format <sys>_induce_defect<N>, return N.
        """
        tempname = self.keywords['name'].lower()
        numstr = tempname[tempname.find("defect")+6:]
        if numstr.find("defect") == -1:
            pass
        else:
            numstr = numstr[numstr.find("defect")+6:] #allow 'defect' to be in <sys>
        defectnum = int(numstr)
        return defectnum

from MAST.ingredients.neb import NEB

class NEBLowMesh(NEB):
    def __init__(self, **kwargs):
        NEB.__init__(self, **kwargs)
    
    def update_children(self):
        """Update to ANOTHER NEB."""
        import os
        for childname in self.keywords['child_dict'].iterkeys():
            myct=1
            while myct <= self.keywords['program_keys']['images']:
                imno = str(myct).zfill(2)
                impath = os.path.join(self.keywords['name'], imno)
                self.forward_parent_structure(impath, childname,"parent_structure_" + '-'.join(self.labels) + '_' + imno)
                myct = myct + 1
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import MASTObj
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.utility import MASTError
from MAST.utility.mastfile import MASTFile
from MAST.utility.metadata import Metadata

import os
import shutil
#import pdb
#TTM

class NEB(BaseIngredient):
    """
        Attributes:
        self.labels <list of str>: list of labels labeling the NEB.
                                    self.labels[0] is the initial state label
                                    self.labels[1] is the final state label
        self.neblines <list of str>: list of NEB atom motion lines for this NEB
    """
    def __init__(self, **kwargs):
        #pdb.set_trace()
        allowed_keys = {
                'name' : (str, str(), 'Name of NEB directory'),
                'program': (str, str(), 'DFT program, e.g. "vasp"'),
                'program_keys': (dict, dict(), 'Program-specific keywords'),
                'child_dict':(dict, dict(), 'Dictionary of child names'),
                'structure': (Structure, None, 'Pymatgen-type structure')
                }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        self.neblines = self.keywords['program_keys']['neblines']
        self.labels = "" 

    def is_complete(self):
        return BaseIngredient.is_complete(self)

    def update_children(self):
        """Update children by forwarding parent structure into image folders.
        """
        myct=1
        impath=""
        for childname in sorted(self.keywords['child_dict'].iterkeys()):
            impath = os.path.join(self.keywords['name'], str(myct).zfill(2))
            self.forward_parent_structure(impath, childname)       
            myct = myct + 1

    def write_files(self):
        """Get the parent structures, sort and match atoms, and interpolate.
            Write images to the appropriate folders.
        """
        self.labels = self.get_my_labels()
        parentstructures = self.get_parent_structures()
        parentimagestructures = self.get_parent_image_structures()
        if len(parentimagestructures) == 0:
            image_structures = self.do_interpolation(parentstructures)
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
   
    def get_my_labels(self):
        """For neb in the format <sys>_neb_<N>-<N>_... return the two
            labels identifying the defect sites.
            This should also be in the metadata.txt file as "neblabel=<N>-<N>"
        """
        myname = self.keywords['name']
        mymeta = Metadata(metafile=os.path.join(myname, "metadata.txt"))
        #tempname = myname.lower()
        #lsplit = tempname.split('_') #last underscore segment
        #nebidx = lsplit.index('neb')
        #lseg = lsplit[nebidx + 1] #Follows after "neb"
        myneblabel = mymeta.search_data("neblabel")
        if myneblabel == "":
            raise MASTError(self.__class__.__name__, 
                "No metadata for tag neblabel.")
        lseg = myneblabel[1]
        llist = lseg.split('-')
        if not(len(llist) == 2):
            raise MASTError(self.__class__.__name__, 
                "Error: number of NEB parents <> 2 in " + myname+":"+str(llist))
            return None
        return llist


    def sort_structure_and_neb_lines(self, mystruc, initorfin):
        """Sort the structure in an approved way:
            1. Edge-coorect the fractional structure so for example 
                -0.05 and 0.95 are the same.
            2. Remove lines which closely match the neb moving lines in the 
            NEB section.
            3. Sort the structure by element and coordinate.
            4. Prepend the NEB moving lines to their element sections.
            5. Return the modified structure.

            mystruc <Structure>   : pymatgen structure
            initorfin <int>       : 0 = initial cfg, 1 = final cfg
        """
        from pymatgen.core.structure_modifier import StructureEditor
        from pymatgen.util.coord_utils import find_in_coord_list
        import MAST.data
        atol = 0.1 # Need fairly large tolerance to account for relaxation.
        sortedstruc = mystruc.get_sorted_structure()
        sortedstruc_base = self.keywords['structure'].get_sorted_structure()
        struct_ed = StructureEditor(sortedstruc)
        struct_ed_base = StructureEditor(sortedstruc_base)
        # We are not actually translating the sites; just get them into unit
        # cell
        #TTM DEBUG REMOVE THIS struct_ed.translate_sites(range(0,len(sortedstruc.sites)),np.zeros(3),True)
        nebidx = list()

        elemstarts = self.get_element_indices(sortedstruc)
        mylabels = self.get_my_labels()
        mylabel = '-'.join(mylabels)
        #print "TTM DEBUG: my neb lines", self.neblines[mylabel]
        for nebline in self.neblines[mylabel]:
            usebase=0
            nebdict = self.parse_neb_line(nebline)
            mycoord = nebdict['coord'][initorfin]
            index = find_in_coord_list(sortedstruc.frac_coords, mycoord, atol)
            print 'TTM DEBUG: index: ', index
            if len(index) == 0: #try the base structure, for vacancies
                usebase=1
                index = find_in_coord_list(sortedstruc_base.frac_coords,
                    mycoord, atol)
            if len(index) == 0:
                raise MASTError(self.__class__.__name__, "No coordinate found matching %s" % mycoord)
            if usebase==0:
                nebidx.append(index[0]) #only take first site?
                mysite = sortedstruc.sites[index[0]]
                myelem = MAST.data.atomic_number[mysite.species_string]
                struct_ed.delete_site(index)
                struct_ed.insert_site(elemstarts[myelem], mysite.specie,
                                        mysite.frac_coords)
            else:
                nebidx.append(index[0]) #only take first site?
                mysite = sortedstruc_base.sites[index[0]]
                myelem = MAST.data.atomic_number[mysite.species_string]
                struct_ed_base.delete_site(index)
                struct_ed_base.insert_site(elemstarts[myelem], mysite.specie,
                                        mysite.frac_coords)

        if not len(nebidx) == len(self.neblines[mylabel]):
            raise MASTError(self.__class__.__name__, "Not all NEB lines found.")
        return struct_ed.modified_structure        

    def get_element_indices(self, sortedstruc):
        """From a sorted structure, get the element indices
            Args:
                sortedstruc <Structure>: pymatgen structure sorted by
                                            electronegativity
            Returns:
                elstart <dict>: element index dictionary, in the format:
                            elstart[element number] = starting index, e.g.
                            elstart[24] = 0
        """
        elstart = dict()
        numlist = sortedstruc.atomic_numbers
        listlen = len(numlist)
        idx = 0
        while idx < listlen:
            atomnum = numlist[idx]
            if not atomnum in elstart.keys():
                elstart[atomnum] = idx
            idx = idx + 1
        return elstart

    def parse_neb_line(self, nebline):
        """Parse an NEB atomic movement line which has the following format:
                    Cr, 0 0 0, 0.5 0.5 0.5 in the input file, and
                    ['Cr','0 0 0','0.5 0.5 0.5'] here as nebline
                    element, init_coord, fin_coord
            Args:
                line <list of str>: NEB atomic movement line
            Returns:
                linedict <dict>: 
                                 ['element'] <int> = element's atomic number
                                 ['coord'][0] <np.array> = initial coordinates
                                 ['coord'][1] <np.array> = final coordinates
        """
        nebdict=dict()
        #print "TTM DEBUG: nebline", nebline
        nebelem = str(nebline[0].strip()).title()
        import MAST.data
        nebdict['element'] = MAST.data.atomic_number[nebelem]
        nebdict['coord'] = dict()
        nebdict['coord'][0] = np.array(nebline[1].split(), dtype='float')
        nebdict['coord'][1] = np.array(nebline[2].split(), dtype='float')
        return nebdict
        
    def get_parent_structures(self):
        """Assume that parents have written files
            named 'parent_structure_<N>'. 
            For VASP these are CONTCAR-type files.
            Returns:
                [struct_init, struct_fin]: pymatgen Structure objects
        """
        header = os.path.join(self.keywords['name'], "parent_structure_")
        pfpath_init = header + self.labels[0]
        pfpath_fin = header + self.labels[1]
        if not os.path.isfile(pfpath_init):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_init)
        if not os.path.isfile(pfpath_fin):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_fin)
        struct_init = BaseIngredient.get_structure_from_file(self, pfpath_init)
        struct_fin = BaseIngredient.get_structure_from_file(self, pfpath_fin)
        sorted_init = self.sort_structure_and_neb_lines(struct_init, 0) 
        sorted_fin = self.sort_structure_and_neb_lines(struct_fin, 1)
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
                sorted_im = self.sort_structure_and_neb_lines(struct_im, 0) 
                imstrs.append(sorted_im)
            imct = imct + 1
        if len(imstrs) > 0 and not (len(imstrs) == numim):
            raise MASTError(self.__class__.__name__, "Incomplete number of forwared images found!")
        return imstrs
    def do_interpolation(self, parentstructures):
        """Do interpolation."""
        if parentstructures == None:
            raise MASTError(self.__class__.__name__,"Bad number of parent paths.")
        struct_init = parentstructures[0]
        struct_fin = parentstructures[1]
        structure_list = struct_init.interpolate(struct_fin, int(self.keywords['program_keys']['images']) + 1 )
        return structure_list

    def place_parent_energy_files(self):
        """Assume parents have written files parent_energy_<N>.
            Copy these files into the 00 and 0N directories.
        """
        header = os.path.join(self.keywords['name'], "parent_energy_")
        pfpath1= header + self.labels[0]
        pfpath2= header + self.labels[1]
        pffile1=MASTFile(pfpath1)
        pffile2=MASTFile(pfpath2)
        pffile1.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self,1)) #MASTFile contains directory locking.
        pffile2.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, 2))
        return

    def run(self, mode='serial', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)

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
                    self.forward_extra_restart_files(newname, childname)
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

from MAST.ingredients.optimize import Optimize

class OptimizeDefect(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
from MAST.ingredients.optimize import Optimize

class OptimizeLowMeshDefect(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
from MAST.ingredients.optimize import Optimize

class OptimizeLowMesh(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import MASTObj
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.utility import MASTError
from MAST.utility import dirutil

import os
import shutil
import subprocess
import time
#TA


class Optimize(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'child_dict': (dict, dict(), 'Dictionary of children'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

        #import inspect
        #print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])
        #print 'Program_keys =', self.keywords['program_keys']

    def is_complete(self):
        return BaseIngredient.is_complete(self)

    def update_children(self):
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname)
            
    def write_files(self):
        self.set_up_program_input()
        self.write_submit_script()
            
    def run(self, mode='serial', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)
import numpy as np
import pymatgen
from MAST.ingredients.phononsingle import PhononSingle
from MAST.ingredients.baseingredient import BaseIngredient
import os
import shutil
from MAST.utility import dirutil

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
        self.combine_dynmats()
        shutil.copy(os.path.join(self.keywords['name'],"DYNMAT_combined"),
            os.path.join(self.keywords['name'],"DYNMAT"))
        self.combine_displacements()
        shutil.copy(os.path.join(self.keywords['name'],"XDATCAR_combined"),
            os.path.join(self.keywords['name'],"XDATCAR"))
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_dynmat(self.keywords['name'], childname)
            #Do NOT forward the CONTCAR, since the ending CONTCAR contains a displacement in it. Forward the POSCAR instead.
            self.forward_parent_initial_structure(self.keywords['name'],childname, "POSCAR_prePHON")

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
    def is_complete(self):
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

    def is_ready_to_run(self):
        """Make sure all subfolders are ready to run."""
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

    def run(self):
        """Run all subfolders."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            BaseIngredient.run(self,"serial")
        self.keywords['name']=myname
        return

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
            self.foward_parent_initial_structure(self.keywords['name'],childname, "POSCAR_prePHON")

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

class PhonParse(Optimize):
    """Phonon parser using PHON
    """
    
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
        print "USING PHON PROGRAM. Cite Dario Alfe."
        self.keywords['program'] = 'phon'

from MAST.ingredients.optimize import Optimize

class SinglePoint(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
from MAST.ingredients.neb import NEB
import os 
from MAST.utility.mastfile import MASTFile
from MAST.ingredients.baseingredient import BaseIngredient

class SortNEBStatic(NEB):
    """This class moves static image folders into a new folder. No calculations
        are run.
    """
    def __init__(self, **kwargs):
        NEB.__init__(self, **kwargs)

    def write_files(self):
        NEB.write_files(self)
        dirname = self.keywords['name']
        dirlist = os.listdir(dirname)
        for myfile in dirlist:
           if 'parent_energy' in myfile:
               imfolder = myfile.split('_')[-1]
               if imfolder in dirlist:
                   cpfile = MASTFile(os.path.join(dirname,myfile))
                   cpfile.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, imfolder))


    def is_complete(self):
        dirname = self.keywords['name']
        dirlist = os.listdir(dirname)
        okay=True
        dirct=0
        for mydir in dirlist:
            if os.path.isdir(os.path.join(dirname,mydir)):
                dirct = dirct + 1
                if len(os.listdir(os.path.join(dirname, mydir))) == 2:
                    pass
                else:
                    okay=False
        if dirct == 0:
            return False
        return okay

    def run(self):
        pass
from MAST.ingredients.optimize import Optimize
from MAST.utility.metadata import Metadata
from MAST.utility import MASTError
import os
class StaticDefectForNEB(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        label = self.get_my_label()
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname,"parent_structure_" + label)
            self.forward_parent_energy(self.keywords['name'], childname, "parent_energy_" + label)

    def get_my_label(self):
        """Return the defect label, without the word "defect".
            Returns:
                label <str>: defect label
        """
        myname = self.keywords['name']
        mymeta = Metadata(metafile=os.path.join(myname, "metadata.txt"))
        mylabel = mymeta.search_data("defect_label")
        if mylabel == "":
            raise MASTError(self.__class__.__name__,
                "No metadata file for tag defect.")
        plabel = mylabel[1]
        if 'defect_' in plabel:
            plabel = plabel.split("defect_")[-1]
        return plabel
from MAST.ingredients.optimize import Optimize

class StaticForPhonon(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        Optimize.update_children(self)
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_extra_restart_files(self.keywords['name'], childname)

from MAST.ingredients.optimize import Optimize

class Static(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
