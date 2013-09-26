from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.cifio import CifParser
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata
import os
import shutil
from MAST.ingredients.checker import BaseChecker
class VaspChecker(BaseChecker):
    """VASP checker functions
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)

    def get_structure_from_file(self, filepath):
        """Get structure from file."""
        return Poscar.from_file(filepath, False).structure

    def get_structure_from_directory(self, dirname):
        """Get structure from directory. Preferentially gets CONTCAR first."""
        cpath = os.path.join(dirname, "CONTCAR")
        ppath = os.path.join(dirname, "POSCAR")
        if os.path.isfile(cpath):
            return Poscar.from_file(cpath, False).structure
        elif os.path.isfile(ppath):
            return Poscar.from_file(ppath, False).structure
        else:
            raise MASTError(self.__class__.__name__, "No valid structure in %s" % dirname)

    def forward_parent_dynmat(self, parentpath, childpath, newname1="DYNMAT", newname2="XDATCAR"):
        """Forward the DYNMAT and the XDATCAR."""
        dirutil.lock_directory(childpath)
        shutil.copy(os.path.join(parentpath, "DYNMAT"),os.path.join(childpath, newname1))
        shutil.copy(os.path.join(parentpath, "XDATCAR"),os.path.join(childpath, newname2))
        dirutil.unlock_directory(childpath)
        return


    def forward_parent_structure(self, parentpath, childpath, newname="POSCAR"):
        """Copy CONTCAR to new POSCAR"""
        dirutil.lock_directory(childpath)
        shutil.copy(os.path.join(parentpath, "CONTCAR"),os.path.join(childpath, newname))
        dirutil.unlock_directory(childpath)
        return

    def forward_parent_initial_structure(self, parentpath, childpath, newname="POSCAR"):
        """Copy POSCAR to new POSCAR (for use in phonons)"""
        dirutil.lock_directory(childpath)
        shutil.copy(os.path.join(parentpath, "POSCAR"),os.path.join(childpath, newname))
        dirutil.unlock_directory(childpath)
        return
    def forward_parent_energy(self, parentpath, childpath, newname="OSZICAR"):
        """Copy OSZICAR"""
        dirutil.lock_directory(childpath)
        shutil.copy(os.path.join(parentpath, "OSZICAR"),os.path.join(childpath, newname))
        dirutil.unlock_directory(childpath)
        return

    def images_complete(self, dirname, numim):
        """Check if all images in a VASP NEB calculation are complete.
            dirname = directory housing /00.../0N+1 files; 
                      only checks directories /01.../0N where N is # images
            numim = number of images
        """
        imct=1
        numim = int(numim)
        while imct <= numim:
            num_str = str(imct).zfill(2)
            impath = os.path.join(dirname, num_str)
            try:
                myoutcar = Outcar(os.path.join(impath, "OUTCAR"))
            except (IOError):
                return False
            if not 'User time (sec)' in myoutcar.run_stats.keys():
                return False
            if myoutcar.run_stats['User time (sec)'] > 0:
                #print "image",imct,"complete"
                pass
            else:
                return False
            imct = imct + 1
        return True


    def is_complete(self, dirname):
        """Check if single VASP non-NEB calculation is complete."""
        try:
            myoutcar = Outcar(os.path.join(dirname, "OUTCAR"))
        except (IOError):
            return False

        #hw 04/19/13
        try:
            if myoutcar.run_stats['User time (sec)'] > 0:
                return True
            else:
                return False
        except KeyError:
            return False

    def is_ready_to_run(self, dirname):
        """Check if vasp is ready to run."""
        notready=0
        if not(os.path.isfile(dirname + "/KPOINTS")):
            notready = notready + 1
        if not(os.path.isfile(dirname + "/POTCAR")):
            notready = notready + 1
        if not(os.path.isfile(dirname + "/INCAR")):
            notready = notready + 1
        else: #INCAR exists. Now check POSCARs.
            myincar = MASTFile(dirname + "/INCAR")
            if myincar.get_line_match("IMAGES") == None: 
                if not(os.path.isfile(dirname + "/POSCAR")):
                    notready = notready + 1
            else:
                subdirs = dirutil.walkdirs(dirname)
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

    def _vasp_poscar_setup(self):
        name = self.keywords['name']
        pospath = os.path.join(name, "POSCAR")
        if os.path.isfile(pospath):
            my_poscar = Poscar.from_file(pospath) 
            #parent should have given a structure
        else: #this is an originating run; mast should give it a structure
            my_poscar = Poscar(self.keywords['structure'])
        if 'mast_coordinates' in self.keywords['program_keys'].keys():
            goodstrucs=[my_poscar.structure.copy()]
            coordstrucs=get_coordinates_only_structure_from_input(keywords)
            newstrucs=graft_coordinates_onto_structure(goodstrucs, 
                coordstrucs)
            my_poscar.structure=newstrucs[0].copy()
        dirutil.lock_directory(name)
        my_poscar.write_file(pospath)
        dirutil.unlock_directory(name)
        return my_poscar

    def _vasp_kpoints_setup(self):
        """Parse mast_kpoints string, which should take the format:
            number, number, number designation
            examples: "3x3x3 M", "1x1x1 G". If no designation is given,
            Monkhorst-Pack is assumed.
        """
        name = self.keywords['name']
        if 'mast_kpoints' in self.keywords['program_keys'].keys():
            kptlist = self.keywords['program_keys']['mast_kpoints']
        else:
            raise MASTError("vasp_checker, _vasp_kpoint_setup","k-point instructions need to be set in ingredients keyword mast_kpoints")
        if len(kptlist) == 3:
            desig = "M"
        else:
            desig = kptlist[3].upper()
        if desig == "M":
            my_kpoints = Kpoints.monkhorst_automatic(kpts=(int(kptlist[0]),int(kptlist[1]),int(kptlist[2])),shift=(0,0,0))
        elif desig == "G":
            my_kpoints = Kpoints.gamma_automatic(kpts=(int(kptlist[0]),int(kptlist[1]),int(kptlist[2])),shift=(0,0,0))
        else:
            raise MASTError("vasp_checker, _vasp_kpoint_setup","kpoint designation " + desig + " not recognized.")
        dirutil.lock_directory(name)
        my_kpoints.write_file(name + "/KPOINTS")
        dirutil.unlock_directory(name)
        return my_kpoints

    def _vasp_potcar_setup(self, my_poscar):
        name = self.keywords['name']

        if 'mast_xc' in self.keywords['program_keys'].keys():
            myxc = self.keywords['program_keys']['mast_xc'].upper() #Uppercase
        else:
            raise MASTError("vasp_checker, _vasp_potcar_setup","Exchange correlation functional needs to be specified in ingredients keyword mast_xc")

        if ('mast_pp_setup' in self.keywords['program_keys'].keys()):
            sites = my_poscar.site_symbols
            setup = self.keywords['program_keys']['mast_pp_setup']
            pp_sites = list()
            for site in sites:
                try:
                    pp_sites.append(setup[site])
                except KeyError:
                    pp_sites.append(site)
            my_potcar = Potcar(symbols=pp_sites, functional=myxc, sym_potcar_map=None)
        else:
            my_potcar = Potcar(symbols=my_poscar.site_symbols, functional=myxc, sym_potcar_map=None)

        dirutil.lock_directory(name)
        my_potcar.write_file(name + "/POTCAR")
        dirutil.unlock_directory(name)
        return my_potcar

    def _vasp_incar_get_non_mast_keywords(self):
        """Get the non-VASP keywords and make a dictionary."""
        incar_dict=dict()
        allowedpath = os.path.join(dirutil.get_mast_install_path(), 'MAST',
                        'ingredients','programkeys','vasp_allowed_keywords.py')
        allowed_list = _vasp_incar_get_allowed_keywords(allowedpath)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                keytry = key.upper()
                if not (keytry in allowed_list):
                    print "Ignoring program key %s for INCAR. To allow this keyword, add it to %s" % (keytry, allowedpath)
                else:
                    if type(value)==str and value.isalpha():
                        incar_dict[keytry]=value.capitalize() #First letter cap
                    else:
                        incar_dict[keytry]=value
        return incar_dict

    def _vasp_is_neb(self):
        """Check if ingredient type is an actual NEB run.
            Returns:
                True if NEB, False otherwise.
        """
        metapath = "%s/metadata.txt" % self.keywords['name']
        if not os.path.isfile(metapath): #we are in some sort of sub-ingredient folder masquerading as a separate ingredient
            if not os.path.isfile("%s/metadata.txt" % os.path.dirname(keywords['name'])): #no metadata can be found
                return False 
            else:
                mymeta = Metadata(metafile="%s/metadata.txt" % os.path.dirname(keywords['name']))
        else:
            mymeta = Metadata(metafile=metapath)
        [ingline,ingval]=mymeta.search_data("neb_label")
        if not ingval == None and (not 'stat' in self.keywords['name']):
            return True
        else:
            return False

    def _vasp_incar_get_allowed_keywords(self, allowedpath):
        """Get allowed vasp keywords.
            Args:
                allowedpath <str>: file path for allowed vasp keywords
        """
        allowed = MASTFile(allowedpath)
        allowed_list=list()
        for entry in allowed.data:
            allowed_list.append(entry.strip())
        return allowed_list

    def _vasp_incar_setup(self, my_potcar, my_poscar):
        """Set up the INCAR, including MAGMOM string, ENCUT, and NELECT."""
        name=self.keywords['name']
        myd = dict()
        myd = self._vasp_incar_get_non_mast_keywords()
        if not self._vasp_is_neb():
            try:
                myd.pop("IMAGES")
            except KeyError:
                pass
        if 'mast_multiplyencut' in self.keywords['program_keys'].keys():
            mymult = float(self.keywords['program_keys']['mast_multiplyencut'])
        else:
            mymult = 1.5
        if 'ENCUT' in myd.keys():
            pass
        else:
            myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(my_potcar)*mymult
        if 'mast_setmagmom' in self.keywords['program_keys'].keys():
            magstr = str(self.keywords['program_keys']['mast_setmagmom'])
            magmomstr=""
            maglist = magstr.split()
            numatoms = sum(my_poscar.natoms)
            if len(maglist) < numatoms:
                magct=0
                while magct < len(maglist):
                    magmomstr = magmomstr + str(my_poscar.natoms[magct]) + "*" + maglist[magct] + " " 
                    magct = magct + 1
            else:
                magmomstr = magstr
            myd['MAGMOM']=magmomstr
        if 'mast_charge' in self.keywords['program_keys'].keys():
            myelectrons = vasp_extensions.get_total_electrons(my_poscar, my_potcar)
            newelectrons=0.0
            try:
                adjustment = float(self.keywords['program_keys']['mast_charge'])
            except (ValueError, TypeError):
                raise MASTError("vasp_checker, vasp_incar_setup","Could not parse adjustment")
            #newelectrons = myelectrons + adjustment
            newelectrons = myelectrons - adjustment
            myd['NELECT']=str(newelectrons)
        my_incar = Incar(myd)
        dirutil.lock_directory(name)
        my_incar.write_file(name + "/INCAR")
        dirutil.unlock_directory(name)
        return my_incar

    def set_up_program_input(self):
        myposcar = self._vasp_poscar_setup()
        self._vasp_kpoints_setup()
        mypotcar = self._vasp_potcar_setup(myposcar)
        self._vasp_incar_setup(mypotcar, myposcar)
        return

    def get_path_to_write_neb_parent_energy(self, myname, myimages, parent):
        if parent == 1:
            return os.path.join(myname, "00", "OSZICAR")
        elif parent == 2:
            return os.path.join(myname, str(int(myimages)+1).zfill(2),"OSZICAR")
        elif len(parent) > 1:
            return os.path.join(myname, parent, "OSZICAR")
        else:
            raise MASTError("vasp_checker, get_path_to_write_neb_parent_energy","Parent not specified correctly.")

    def set_up_neb_folders(self, myname, image_structures, keywords):
        imct=0
        if 'mast_coordinates' in keywords['program_keys'].keys():
            goodstrucs=list()
            for imstruc in image_structures[1:-1]:
                goodstrucs.append(imstruc.copy())
                coordstrucs=get_coordinates_only_structure_from_input(keywords)
                newstrucs=graft_coordinates_onto_structure(goodstrucs, 
                    coordstrucs)
        while imct < len(image_structures):
            imposcar = Poscar(image_structures[imct])
            num_str = str(imct).zfill(2)
            impath = os.path.join(myname, num_str)
            impospath = os.path.join(myname, "POSCAR_" + num_str)
            if 'mast_coordinates' in keywords['program_keys'].keys():
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
                print "Directory at", impath, "already exists."
                return None
            dirutil.lock_directory(impath)
            imposcar.write_file(os.path.join(impath, "POSCAR"))
            dirutil.unlock_directory(impath)
            imct = imct + 1
        return
        

    def set_up_program_input_neb(self, keywords, image_structures):
        set_up_neb_folders(keywords['name'], image_structures, keywords)
        _vasp_kpoints_setup(keywords)
        mypotcar = _vasp_potcar_setup(keywords, Poscar(image_structures[0]))
        _vasp_incar_setup(keywords, mypotcar, Poscar(image_structures[0]))
        return
    def forward_extra_restart_files(self, parentpath, childpath):
        """Forward extra restart files: softlink to WAVECAR and CHGCAR."""
        dirutil.lock_directory(childpath)
        import subprocess
        #print "cwd: ", os.getcwd()
        curpath = os.getcwd()
        os.chdir(childpath)
        #print parentpath
        #print childpath
        mylink=subprocess.Popen("ln -s %s/WAVECAR WAVECAR" % parentpath, shell=True)
        mylink.wait()
        mylink2=subprocess.Popen("ln -s %s/CHGCAR CHGCAR" % parentpath, shell=True)
        mylink2.wait()
        os.chdir(curpath)
        dirutil.unlock_directory(childpath)
    def add_selective_dynamics_to_structure(keywords, sdarray):
        name = keywords['name']
        pname = os.path.join(name,"POSCAR")
        phposcar = Poscar.from_file(pname)
        phposcar.selective_dynamics = sdarray
        dirutil.lock_directory(name)
        os.rename(pname, pname + "_no_sd")
        phposcar.write_file(pname)
        dirutil.unlock_directory(name)
        return

    def get_vasp_energy(abspath):
        return Vasprun('%s/vasprun.xml' % abspath).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

    def get_coordinates_only_structure_from_input(self):
        """Get coordinates-only structures from mast_coordinates
            ingredient keyword
            Args:
                keywords <dict>: ingredient keywords
            Returns:
                coordstrucs <list>: list of Structure objects
        """
        coordposlist=self.keywords['program_keys']['mast_coordinates']
        coordstrucs=list()
        coordstruc=None
        for coordpositem in coordposlist:
            if ('poscar' in os.path.basename(coordpositem).lower()):
                coordstruc = Poscar.from_file(coordpositem).structure
            elif ('cif' in os.path.basename(coordpositem).lower()):
                coordstruc = CifParser(coordpositem).get_structures()[0]
            else:
                error = 'Cannot build structure from file %s' % coordpositem
                raise MASTError("vasp_checker,get_coordinates_only_structure_from_input", error)
            coordstrucs.append(coordstruc)
        return coordstrucs

    def graft_coordinates_onto_structure(self, goodstrucs,coordstrucs):
        """Graft coordinates from mast_coordinates Structure objects
            onto the appropriate structure
            Args:
                goodstrucs <list>: list of Structure objects with
                    the correct lattice parameters and elements
                coordstrucs <list>: list of Structure objects with
                    the coordinates for grafting
            Returns:
                modstrucs <list>: list of modified Structure objects
        """
        sidx=0
        slen=len(goodstrucs)
        while sidx < slen:
            lengoodsites=len(goodstrucs[sidx].sites)
            lencoordsites=len(coordstrucs[sidx].sites)
            if not (lengoodsites == lencoordsites):
                raise MASTError("vasp_checker,graft_coordinates_onto_structure", "Original and coordinate structures do not have the same amount of sites.")
            cct=0
            newsites=list()
            mylattice=goodstrucs[sidx].lattice
            while cct < lengoodsites:
                newcoords=coordstrucs[sidx].sites[cct].frac_coords
                oldspecie=goodstrucs[sidx].sites[cct].specie
                newsite=PeriodicSite(oldspecie, newcoords, mylattice)
                newsites.append(newsite)
                cct=cct+1
            goodstrucs[sidx].remove_sites(range(0,lengoodsites))
            for cct in range(0, lengoodsites):
                goodstrucs[sidx].append(newsites[cct].specie,
                    newsites[cct].frac_coords)
            sidx = sidx + 1
        return goodstrucs

    def get_max_enmax_from_potcar(mypotcar):
        """Get maximum enmax value (float) from Potcar (combined list)"""
        enmax_list=list()
        potcarct=0
        onepotcar=None
        while potcarct < len(mypotcar):
            onepotcar = mypotcar[potcarct] #A PotcarSingle object
            enmax_list.append(onepotcar.enmax)
            potcarct = potcarct + 1
        return max(enmax_list)

    def make_one_unfrozen_atom_poscar(myposcar, natom):
        """Use selective dynamics to make a poscar with one unfrozen atom.
            myposcar = Poscar
            natom = the number of the atom to unfreeze
            Returns: Poscar (use write_file function on it).
        """
        mysd=np.zeros([sum(myposcar.natoms),3],bool)
        mysd[natom-1][0]=True #indexing starts at 0
        mysd[natom-1][1]=True
        mysd[natom-1][2]=True
        myposcar.selective_dynamics = mysd
        return myposcar

    def make_one_unfrozen_direction_poscar(myposcar, natom, ndir):
        """Use selective dynamics to make a poscar with one unfrozen atom.
            myposcar = Poscar
            natom = the number of the atom to unfreeze
            ndir = the direction to freeze (0, 1, 2 for x, y, z)
            Returns: Poscar (use write_file function on it).
        """
        mysd=np.zeros([sum(myposcar.natoms),3],bool)
        mysd[natom-1][ndir]=True #indexing starts at 0
        myposcar.selective_dynamics = mysd
        return myposcar

        """Recombine DYNMAT files from one unfrozen atom in folder."""
        natoms = sum(myposcar.natoms)
        myhess=np.zeros([natoms*3, natoms*3])
        #arrange as x1 y1 z1 x2 y2 z2 x3 y3 z3...
        dirlist = os.listdir(mydir)
        dynlist=[]
#    for entry in dirlist:

    def read_my_dynmat(mydir, fname="DYNMAT"):
        """Read a DYNMAT file.
            Returns:
                dyndict <dict>: dictionary structured like this:
                    ['numspec'] = <int> number of species
                    ['numatoms'] = <int> number of atoms
                    ['numdisp'] = <int> number of displacements
                    ['massline'] = <str> masses line
                    ['atoms'][atom <int>][disp <int>]['displine'] =
                        displacement
                        vector (part of first line in dynmat block,
                        e.g. "0.01 0 0")
                    ['atoms'][atom <int>][disp <int>]['dynmat'] =
                            <list>
                            list of dynmat lines for this atom
                            and this displacement
        """
        mydyn = MASTFile(os.path.join(mydir, fname))
        mydata = list(mydyn.data) #whole new copy
        dyndict=dict()
        firstline = mydata.pop(0) #pop first value
        firstspl = firstline.strip().split()
        dyndict['numspec']=int(firstspl[0])
        numatoms = int(firstspl[1])
        dyndict['numatoms']=numatoms
        dyndict['numdisp']=int(firstspl[2])
        dyndict['massline']=mydata.pop(0)
        dyndict['atoms']=dict()
        atom=0
        disp=0
        thirdline=""
        while len(mydata) > 0:
            thirdline = mydata.pop(0)
            thirdspl = thirdline.strip().split()
            atom=int(thirdspl[0])
            disp=int(thirdspl[1])
            displine=' '.join(thirdspl[2:])
            if not atom in dyndict['atoms'].keys():
                dyndict['atoms'][atom]=dict()
            dyndict['atoms'][atom][disp]=dict()
            dyndict['atoms'][atom][disp]['displine'] = displine
            dyndict['atoms'][atom][disp]['dynmat']=list()
            for act in range(0, numatoms):
                dyndict['atoms'][atom][disp]['dynmat'].append(mydata.pop(0))
        return dyndict

    def write_my_dynmat(mydir, dyndict, fname="DYNMAT"):
        """Write a DYNMAT file.
            Args:
                mydir <str>: Directory in which to write
                dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
                fname <str>: filename (default DYNMAT)
        """
        dynwrite=MASTFile()
        dynwrite.data=list()
        firstline=str(dyndict['numspec']) + " " + str(dyndict['numatoms']) + " " + str(dyndict['numdisp']) + "\n"
        dynwrite.data.append(firstline)
        dynwrite.data.append(dyndict['massline'])
        atomlist=dyndict['atoms'].keys()
        atomlist.sort()
        for atom in atomlist:
            displist = dyndict['atoms'][atom].keys()
            displist.sort()
            for disp in displist:
                thirdline = str(atom) + " " + str(disp) + " " + dyndict['atoms'][atom][disp]['displine'] + "\n"
                dynwrite.data.append(thirdline)
                for line in dyndict['atoms'][atom][disp]['dynmat']:
                    dynwrite.data.append(line)
        dynwrite.to_file(os.path.join(mydir, fname))


    def write_my_dynmat_without_disp_or_mass(mydir, dyndict, fname="DYNMAT"):
        """Write a DYNMAT file without the displacement indicators 1, 2, 3
            and without the masses line, and with first line having only
            the total number of displacements, for PHON.
            Args:
                mydir <str>: Directory in which to write
                dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
                fname <str>: filename (default DYNMAT)
        """
        dynwrite=MASTFile()
        dynwrite.data=list()
        firstline=str(dyndict['numdisp']) + "\n"
        dynwrite.data.append(firstline)
        atomlist=dyndict['atoms'].keys()
        atomlist.sort()
        for atom in atomlist:
            displist = dyndict['atoms'][atom].keys()
            displist.sort()
            for disp in displist:
                thirdline = str(atom) + " " + dyndict['atoms'][atom][disp]['displine'] + "\n"
                dynwrite.data.append(thirdline)
                for line in dyndict['atoms'][atom][disp]['dynmat']:
                    dynwrite.data.append(line)
        dynwrite.to_file(os.path.join(mydir, fname))

    def read_my_xdatcar(mydir, fname="XDATCAR"):
        """Read an XDATCAR file.
            Returns:
                xdatdict <dict>: Dictionary of configurations
                    xdatdict['descline'] = <str> description line
                    xdatdict['specline'] = <str> species line
                    xdatdict['numline'] = <str> numbers lines
                    xdatdict['type'] = <str> Direct or not
                    xdatdict['numatoms'] = <int> number of atoms
                    xdatdict['configs'][config <int>] = <list>
                            list of lines in this configuration
                    xdatdict['scale'] = <str> scale
                    xdatdict['latta'] = <str> lattice vector a
                    xdatdict['lattb'] = <str> lattice vector b
                    xdatdict['lattc'] = <str> lattice vector c
        """
        myxdat = MASTFile(os.path.join(mydir, fname))
        mydata = list(myxdat.data) #whole new copy
        xdatdict=dict()
        xdatdict['descline'] = mydata.pop(0) #pop first value
        tryspec = mydata.pop(0)
        tryspecstrip = tryspec.strip()
        tryspecstrip = tryspecstrip.replace(" ","x")
        if not tryspecstrip.isalpha(): # version 5.2.11 and on
            xdatdict['scale'] = tryspec
            xdatdict['latta'] = mydata.pop(0)
            xdatdict['lattb'] = mydata.pop(0)
            xdatdict['lattc'] = mydata.pop(0)
            xdatdict['specline'] = mydata.pop(0)
        else:
            xdatdict['scale'] = ""
            xdatdict['latta'] = ""
            xdatdict['lattb'] = ""
            xdatdict['lattc'] = ""
            xdatdict['specline'] = tryspec
        numline = mydata.pop(0)
        xdatdict['numline'] = numline
        numatoms = sum(map(int, numline.strip().split()))
        xdatdict['numatoms'] = numatoms
        if not tryspecstrip.isalpha(): #version 5.2.11 and on
            xdatdict['type'] = ""
        else:
            xdatdict['type'] = mydata.pop(0)
        xdatdict['configs']=dict()
        kfgct=1
        while len(mydata) > 0:
            mydata.pop(0) # Konfig line, or a blank line
            xdatdict['configs'][kfgct] = list()
            for act in range(0,numatoms):
                if len(mydata) > 0:
                    xdatdict['configs'][kfgct].append(mydata.pop(0))
            kfgct=kfgct+1
        return xdatdict

    def write_my_xdatcar(mydir, xdatdict, fname="XDATCAR"):
        """Write an XDATCAR file.
            Args:
                mydir <str>: Directory in which to write
                xdatdict <dict>: Dictionary of XDATCAR (see read_my_xdatcar)
                fname <str>: filename (default DYNMAT)
        """
        xdatwrite=MASTFile()
        xdatwrite.data=list()
        xdatwrite.data.append(xdatdict['descline'])
        if not (xdatdict['scale'] == ""):
            xdatwrite.data.append(xdatdict['scale'])
            xdatwrite.data.append(xdatdict['latta'])
            xdatwrite.data.append(xdatdict['lattb'])
            xdatwrite.data.append(xdatdict['lattc'])
        xdatwrite.data.append(xdatdict['specline'])
        xdatwrite.data.append(xdatdict['numline'])
        if not (xdatdict['type'] == ""):
            xdatwrite.data.append(xdatdict['type'])
        configlist = xdatdict['configs'].keys()
        configlist.sort()
        for cfg in configlist:
            xdatwrite.data.append("Konfig=%1i\n" % cfg)
            xdatwrite.data.extend(xdatdict['configs'][cfg])
        xdatwrite.to_file(os.path.join(mydir, fname))

    def combine_dynmats(mydir):
        """Combine DYNMATs into one file.
            Args:
                mydir <str>: top directory for DYNMAT files
        """
        dynmatlist = walkfiles(mydir, 2, 5, "*DYNMAT*") #start one level below
        if len(dynmatlist) == 0:
            raise MASTError("pmgextend combine_dynmats", "No DYNMATs found under " + mydir)
        totnumdisp=0
        largedyn=dict()
        largedyn['atoms'] = dict()
        for onedynmat in dynmatlist:
            dyndir = os.path.dirname(onedynmat)
            onedyn = read_my_dynmat(dyndir)
            totnumdisp = totnumdisp + onedyn['numdisp']
            for atom in onedyn['atoms'].keys():
                if not atom in largedyn['atoms'].keys():
                    largedyn['atoms'][atom]=dict()
                    mydisp=1 #start at 1
                for disp in onedyn['atoms'][atom].keys():
                    if disp in largedyn['atoms'][atom].keys():
                        mydisp=mydisp + 1 #increment
                    largedyn['atoms'][atom][mydisp]=dict()
                    largedyn['atoms'][atom][mydisp]['displine'] = str(onedyn['atoms'][atom][disp]['displine'])
                    largedyn['atoms'][atom][mydisp]['dynmat']=list(onedyn['atoms'][atom][disp]['dynmat'])
        largedyn['numspec'] = onedyn['numspec'] #should all be the same
        largedyn['numatoms'] = onedyn['numatoms']
        largedyn['massline'] = onedyn['massline']
        largedyn['numdisp'] = totnumdisp
        write_my_dynmat(mydir, largedyn, "DYNMAT_combined")

    def combine_displacements(mydir):
        """Combine displacements (here XDATCARs) into one file.
            Args:
                mydir <str>: top directory for DYNMAT files
        """
        largexdat=dict()
        xdatlist = walkfiles(mydir, 2, 5, "*XDATCAR*") #start one level below
        if len(xdatlist) == 0:
            raise MASTError("pmgextend combine_displacements", "No XDATCARs found under " + mydir)
        configs=list()
        kfgct=1 # skip config 1 until the end
        largexdat['configs']=dict()
        xdatlist.sort() #get them in order
        for onexdatmat in xdatlist:
            xdatdir = os.path.dirname(onexdatmat)
            onexdat = read_my_xdatcar(xdatdir)
            for kfg in onexdat['configs'].keys():
                if kfg == 1: #skip config 1 until the end
                    pass
                else:
                    kfgct = kfgct + 1 #start at 2
                    largexdat['configs'][kfgct] = onexdat['configs'][kfg]
        largexdat['configs'][1] = onexdat['configs'][1] #get one of the first configs (the first should be the same for all the XDATCARs)
        largexdat['descline'] = onexdat['descline'] #should all be the same
        largexdat['specline'] = onexdat['specline']
        largexdat['scale'] = onexdat['scale']
        largexdat['latta'] = onexdat['latta']
        largexdat['lattb'] = onexdat['lattb']
        largexdat['lattc'] = onexdat['lattc']
        largexdat['numline'] = onexdat['numline']
        largexdat['numatoms'] = onexdat['numatoms']
        largexdat['type'] = onexdat['type']
        write_my_xdatcar(mydir, largexdat, "XDATCAR_combined")
    def make_hessian(myposcar, mydir):
        """Combine DYNMATs into one hessian and solve for frequencies.
            myposcar = Poscar
            mydir = top directory for DYNMAT files
        """
        natoms = sum(myposcar.natoms)
        myhess=np.zeros([natoms*3, natoms*3])
        #arrange as x1, y1, z1, x2, y2, z2, etc.
        dynmatlist = walkfiles(mydir, 1, 5, "*DYNMAT*")
        if len(dynmatlist) == 0:
            raise MASTError("pmgextend combine_dynmats", "No DYNMATs found under " + mydir)
        opendyn=""
        print "DYNMATLIST:"
        print dynmatlist
        datoms=0
        dmats=0
        mct=0
        for onedynmat in dynmatlist:
            dynlines=[]
            opendyn = open(onedynmat,'rb')
            dynlines=opendyn.readlines()
            opendyn.close()
            datoms = int(dynlines[0].split()[1])
            dmats = int(dynlines[0].split()[2])
            mycount=2 #starting line
            while mycount < len(dynlines)-datoms:
                littlemat=[]
                topatom = int(dynlines[mycount].split()[0])
                whichdir = int(dynlines[mycount].split()[1])
                littlemat = dynlines[mycount+1:mycount+datoms+1]
                print littlemat
                act = 0
                while act < datoms:
                    dactx=float(littlemat[act].split()[0])
                    dacty=float(littlemat[act].split()[1])
                    dactz=float(littlemat[act].split()[2])
                    colidx = (topatom-1)*3 + (whichdir-1)
                    #so 2  3  on first line means atom 2's z direction
                    #then with atom 1x, 1y, 1z; 2x, 2y, 2z, etc.
                    myhess[colidx][act*3+0]=dactx
                    myhess[colidx][act*3+1]=dacty
                    myhess[colidx][act*3+2]=dactz
                    act = act + 1
                mycount = mycount + datoms + 1
                print mycount
        print "UNALTERED HESSIAN:"
        print(myhess)
        #create mass matrix
        masses=dynlines[1].split()
        print "MASSES:", masses
        massarr=np.zeros([datoms*3,1])
        act=0
        print myposcar.natoms
        nspec=len(myposcar.natoms)
        totatoms=0
        while act < datoms:
            mymass=0
            nct=0
            totatoms=0
            while (mymass==0) and nct < nspec:
                totatoms = totatoms + myposcar.natoms[nct]
                if act < totatoms:
                    mymass = float(masses[nct])
                nct = nct + 1
            print mymass
            massarr[act*3+0][0]=mymass
            massarr[act*3+1][0]=mymass
            massarr[act*3+2][0]=mymass
            act = act + 1
        massmat = massarr*np.transpose(massarr)
        print "MASS MAT:"
        print massmat
        print "STEP:"
        step = float(dynlines[2].split()[2])
        print step
        print "HESSIAN * -1 / step / sqrt(mass1*mass2)"
        normhess=np.zeros([natoms*3,natoms*3])
        cidx=0
        while cidx < natoms*3:
            ridx=0
            while ridx < natoms*3:
                normhess[ridx][cidx]=-1*myhess[ridx][cidx]/step/np.sqrt(massmat[ridx][cidx])
                ridx = ridx + 1
            cidx = cidx + 1
        print normhess
        print "EIGENVALUES:"
        myeig = np.linalg.eig(normhess)[0]
        print myeig
        print "SQRT of EIGENVALUES in sqrt(eV/AMU)/Angstrom/2pi:"
        myfreq = np.sqrt(myeig)
        print myfreq
        print "SQRT OF EIGENVALUES in THz:"
        myfreqThz = myfreq*15.633302
        print myfreqThz
        myfreqThzsorted = myfreqThz
        myfreqThzsorted.sort()
        print myfreqThzsorted
        return myfreqThzsorted

    def get_total_electrons(myposcar, mypotcar):
        """Get the total number of considered electrons in the system."""
        atomlist = myposcar.natoms
        zvallist = get_zval_list(mypotcar)
        totzval = 0.0
        atomct = 0
        if not (len(zvallist) == len(atomlist)):
            raise MASTError("pmgextend, get_total_electrons",
                "Number of species and number of POTCARs do not match.")
        while atomct < len(atomlist):
            totzval = totzval + (atomlist[atomct] * zvallist[atomct])
            atomct = atomct + 1
        return totzval

    def get_zval_list(mypotcar):
        """Get zvals from POTCAR"""
        zval_list=list()
        potcarct=0
        onepotcar=None
        while potcarct < len(mypotcar):
            onepotcar = mypotcar[potcarct] #A PotcarSingle object
            zval_list.append(onepotcar.zval)
            potcarct = potcarct + 1
        return zval_list
    def get_e0_energy(mydir):
        """Get last E0 energy from OSZICAR.
            Args:
                mydir <str>: Directory in which to look.
            Returns:
                <float>: last E0 energy from OSZICAR
        """
        fullpath=os.path.join(mydir, "OSZICAR")
        if not os.path.isfile(fullpath):
            raise MASTError("vasp_checker, get_e0_energy", "No OSZICAR file at %s" % mydir)
        myosz = MASTFile(fullpath)
        mye0 = myosz.get_segment_from_last_line_match("E0", "E0=","d E =")
        return float(mye0)

