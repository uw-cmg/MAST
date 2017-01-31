##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
from pymatgen.io.vasp import Poscar, Outcar, Potcar, Incar, Kpoints, Vasprun
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
from MAST.ingredients.pmgextend.atom_index import AtomIndex
import os
import shutil
import logging
import subprocess
import numpy as np
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
        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])

    def get_structure_from_file(self, myfilepath=""):
        """Get the structure from a specified file path.
            For VASP, this is a POSCAR-type file.
            Args:
                myfilepath <str>: File path for structure.
        """
        return Poscar.from_file(myfilepath).structure

    def get_initial_structure_from_directory(self,mydir=""):
        """Get the initial structure.
            For VASP, this is the POSCAR file.
            Args:
                mydir <str>: Directory. If no directory is given,
                    use current ingredient directory given as
                    keyword 'name'
        """
        if mydir == "":
            mydir = self.keywords['name']
        return Poscar.from_file("%s/POSCAR" % mydir).structure
    
    def get_final_structure_from_directory(self, mydir=""):
        """Get the final structure.
            For VASP, this is the CONTCAR file.
            Args:
                mydir <str>: Directory. If no directory is given,
                    use current ingredient directory given
                    as keyword 'name'
        """
        if mydir == "":
            mydir = self.keywords['name']
        return Poscar.from_file("%s/CONTCAR" % mydir).structure
    
    def forward_final_structure_file(self, childpath, newname="POSCAR"):
        """Forward the final structure.
            For VASP, this is the CONTCAR.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'POSCAR')
        """
        proceed=False
        workdir=os.path.dirname(self.keywords['name'])
        sdir=os.path.join(workdir,"structure_index_files")
        if os.path.exists(sdir):
            proceed=True
        if not proceed:
            return self.copy_a_file(childpath, "CONTCAR", newname)
        childmeta=Metadata(metafile="%s/metadata.txt" % childpath)
        child_program=childmeta.read_data("program")
        if not "vasp" in child_program: #madelung utility or another folder
            return self.copy_a_file(childpath, "CONTCAR", newname)
        child_scaling_label=childmeta.read_data("scaling_label")
        child_defect_label=childmeta.read_data("defect_label")
        child_neb_label=childmeta.read_data("neb_label")
        child_phonon_label=childmeta.read_data("phonon_label")
        if child_scaling_label == None:
            child_scaling_label = ""
        if child_defect_label == None:
            child_defect_label = ""
        if child_neb_label == None:
            child_neb_label = ""
        if child_phonon_label == None:
            child_phonon_label = ""
        parentmeta=Metadata(metafile="%s/metadata.txt" % self.keywords['name'])
        parent_defect_label=parentmeta.read_data("defect_label")
        parent_neb_label=parentmeta.read_data("neb_label")
        if parent_defect_label == None:
            parent_defect_label = ""
        if parent_neb_label == None:
            parent_neb_label = ""
        if (not (child_neb_label == "")) and (not (parent_defect_label == "")):
            child_defect_label = parent_defect_label
        if (not (child_phonon_label == "")):
            if (not (parent_defect_label == "")):
                child_defect_label = parent_defect_label
            if (not (parent_neb_label == "")):
                child_neb_label = parent_neb_label
                child_defect_label = parent_neb_label.split('-')[0].strip() # always define phonons from first endpoint
        #get child manifest
        childmanifest="manifest_%s_%s_%s" % (child_scaling_label, child_defect_label, child_neb_label)
        #build structure from atom indices using parent name_frac_coords
        ing_label=os.path.basename(self.keywords['name'])
        childmeta.write_data("parent",ing_label)
        mystr=Poscar.from_file("%s/CONTCAR" % self.keywords['name']).structure
        myatomindex=AtomIndex(structure_index_directory=sdir)
        if "inducescaling" in childpath: #initial scaled coords have no parent
            ing_label="original"
        newstr=myatomindex.graft_new_coordinates_from_manifest(mystr, childmanifest,ing_label)
        newposcar=Poscar(newstr)
        self.write_poscar_with_zero_velocities(newposcar, os.path.join(childpath, newname))
        return

    def update_atom_index_for_complete(self):
        """Update atom index files with positions for the 
            completed ingredient.
        """
        proceed=False
        mydir = self.keywords['name']
        workdir=os.path.dirname(mydir)
        sdir=os.path.join(workdir,"structure_index_files")
        if os.path.exists(sdir):
            proceed=True
        if not proceed:
            self.logger.warning("Called update atom index for ingredient %s, but no atom indices exist." % self.keywords['name'])
            return
        mymeta=Metadata(metafile="%s/metadata.txt" % mydir)
        phonon_label = mymeta.read_data("phonon_label")
        if not (phonon_label == None): #Don't update atom indices for phonons
            self.logger.debug("Skipping atom index update for phonon calculation %s" % self.keywords['name'])
            return 
        scaling_label=mymeta.read_data("scaling_label")
        defect_label=mymeta.read_data("defect_label")
        neb_label=mymeta.read_data("neb_label")
        if scaling_label == None:
            scaling_label = ""
        if defect_label == None:
            defect_label = ""
        if neb_label == None:
            neb_label = ""
        if neb_label == "":
            mystr=Poscar.from_file("%s/CONTCAR" % self.keywords["name"]).structure
            ing_label=os.path.basename(mydir)
            manname="manifest_%s_%s_%s" % (scaling_label, defect_label, neb_label)
            myatomindex=AtomIndex(structure_index_directory=sdir)
            myatomindex.update_atom_indices_from_structure(mystr, ing_label, manname)
        else:
            nebdirs = dirutil.immediate_subdirs(self.keywords["name"])
            for mysubdir in nebdirs:
                if os.path.isfile("%s/%s/CONTCAR" % (self.keywords["name"],mysubdir)):
                    mystr=Poscar.from_file("%s/%s/CONTCAR" % (self.keywords["name"], mysubdir)).structure
                    ing_label=os.path.join(os.path.basename(mydir), mysubdir)
                    for defect_label in neb_label.split("-"):
                        manname="manifest_%s_%s_%s" % (scaling_label, defect_label, neb_label)
                        myatomindex=AtomIndex(structure_index_directory=sdir)
                        myatomindex.update_atom_indices_from_structure(mystr, ing_label, manname)
                
        return
    def write_poscar_with_zero_velocities(self, pmg_Poscar, fpath):
        """Write a POSCAR-type file but set zero velocities
            Args:
                pmg_Poscar: pymatgen Poscar objecdt
                fpath: filename
        """
        mystruc = pmg_Poscar.structure
        myvels = np.zeros((mystruc.num_sites,3),'float')
        pmg_Poscar.velocities=myvels
        pmg_Poscar.write_file(fpath)
        return
    def forward_initial_structure_file(self, childpath, newname="POSCAR"):
        """Forward the initial structure.
            For VASP, this is the POSCAR. This function is
            used after phonon calculations, where the CONTCAR
            contains the last displacement. To forward to PHON,
            the POSCAR (without displacements) should be used.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'POSCAR')
        """
        return self.copy_a_file(childpath, "POSCAR", newname)

    def forward_dynamical_matrix_file(self, childpath, newname="DYNMAT"):
        """Forward the dynamical matrix.
            For VASP, this is the DYNMAT file.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'DYNMAT')
        """
        return self.copy_a_file(childpath, "DYNMAT", newname)

    def forward_displacement_file(self, childpath, newname="XDATCAR"):
        """Forward displacement information.
            For VASP, this is the XDATCAR file.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'XDATCAR')
        """
        return self.copy_a_file(childpath, "XDATCAR", newname)

    def forward_energy_file(self, childpath, newname="OSZICAR"):
        """Forward the energy file.
            For VASP, this is the OSZICAR file.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'OSZICAR')
        """
        return self.copy_a_file(childpath, "OSZICAR", newname)

    def is_frozen(self):
        """Check if single VASP non-NEB calculation is frozen.
        """
        return BaseChecker.is_frozen(self, "OUTCAR")

    def is_complete(self):
        """Check if single VASP non-NEB calculation is complete.
        """
        usertime=False
        reachedaccuracy=False
        opath = os.path.join(self.keywords['name'],"OUTCAR")
        if not os.path.isfile(opath):
            self.logger.info("No OUTCAR at %s; not complete." % opath)
            return False

        myoutcar = Outcar(opath)
        
        #hw 04/19/13
        try:
            if myoutcar.run_stats['User time (sec)'] > 0:
                usertime=True
            else:
                usertime=False
        except KeyError:
            usertime=False
        
        reachgrep=subprocess.Popen('grep "reached required accuracy" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            reachedaccuracy=False
        else:
            reachedaccuracy=True

        isstatic=False
        isphonon=False
        isMD=False
        if 'ibrion' in self.keywords['program_keys'].keys():
            ibrionval = str(self.keywords['program_keys']['ibrion'])
            if ibrionval == "-1":
                self.logger.info("ingredient %s IBRION -1 is static" % self.keywords['name'])
                isstatic = True
            elif ibrionval == "0":
                self.logger.info("ingredient %s IBRION 0 is MD" % self.keywords['name'])
                isMD = True
            elif ibrionval in ["5","6","7","8"]:
                self.logger.info("ingredient %s IBRION %s is phonon" % (self.keywords['name'],ibrionval))
                isphonon = True
        if 'nsw' in self.keywords['program_keys'].keys():
            nswval = str(self.keywords['program_keys']['nsw'])
            if nswval in ["0","1"]:
                self.logger.info("ingredient %s NSW %s is static" % (self.keywords['name'],nswval))
                isstatic = True
        else:
            self.logger.info("ingredient %s no NSW given is static" % self.keywords['name'])
            isstatic = True #No NSW specified, VASP defaults to 0
        
        #For static runs, make an additional check for just electronic convergence.
        if isstatic:
            reachgrep=subprocess.Popen('grep "EDIFF is reached" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            reachrpt=reachgrep.communicate()[0]
            reachgrep.wait()
            if reachrpt=='':
                reachedaccuracy=False
            else:
                reachedaccuracy=True


        if (isMD) or (isphonon):
            if not usertime:
                self.logger.warning("OUTCAR at %s shows no user time." % opath)
                return False
            else:
                self.logger.info("OUTCAR at %s shows user time." % opath)
                return True

        if isstatic:
            if reachedaccuracy:
                if usertime:
                    self.logger.info("OUTCAR at %s shows EDIFF reached for static run and user time; complete." % opath)
                    return True
                else:
                    self.logger.warning("OUTCAR at %s shows EDIFF reached for static run, but no user time." % opath)
                    return False
            else:
                if usertime:
                    self.logger.error("OUTCAR at %s does not show EDIFF reached for static run, but shows complete." % opath)
                    return False
                else:
                    self.logger.info("OUTCAR at %s does not show EDIFF reached for static run or user time; still incomplete." % opath)
                    return False

        else:
            if reachedaccuracy:
                if usertime:
                    self.logger.info("OUTCAR at %s shows reached required accuracy and user time; complete." % opath)
                    return True
                else:
                    self.logger.warning("OUTCAR at %s shows reached required accuracy, but no user time." % opath)
                    return False
            else:
                if usertime:
                    self.logger.error("OUTCAR at %s does not show reached required accuracy, but shows complete." % opath)
                    return False
                else:
                    self.logger.info("OUTCAR at %s does not show reached required accuracy or user time; still incomplete." % opath)
                    return False

    def is_ready_to_run(self):
        """Check if single VASP non-NEB ingredient is 
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
        if not(os.path.isfile(dirname + "/POSCAR")):
            notready = notready + 1
        if not(os.path.isfile(dirname + "/submit.sh")):
            notready = notready + 1
        if notready > 0:
            return False
        else:
            return True

    def _vasp_poscar_setup(self):
        """Set up the POSCAR file for a single VASP run.
        """
        name = self.keywords['name']
        pospath = os.path.join(name, "POSCAR")
        if os.path.isfile(pospath):
            my_poscar = Poscar.from_file(pospath) 
            #parent should have given a structure
        else: #this is an originating run; mast should give it a structure
            my_poscar = Poscar(self.keywords['structure'])
            workdir=os.path.dirname(name)
            sdir=os.path.join(workdir,"structure_index_files")
            if os.path.exists(sdir):
                mystr=my_poscar.structure
                manname="manifest___"
                myatomindex=AtomIndex(structure_index_directory=sdir)
                newstr=myatomindex.graft_new_coordinates_from_manifest(mystr, manname, "")
                self.logger.info("Getting original coordinates from manifest.")
                new_pos=Poscar(newstr)
                my_poscar=new_pos
            self.logger.info("No POSCAR found from a parent; base structure used for %s" % self.keywords['name'])
        if 'mast_coordinates' in self.keywords['program_keys'].keys():
            sxtend = StructureExtensions(struc_work1=my_poscar.structure, name=self.keywords['name'])
            coordstrucs=self.get_coordinates_only_structure_from_input()
            newstruc = sxtend.graft_coordinates_onto_structure(coordstrucs[0])
            my_poscar.structure=newstruc.copy()
        dirutil.lock_directory(name)
        self.write_poscar_with_zero_velocities(my_poscar, pospath)
        dirutil.unlock_directory(name)
        return my_poscar

    def _vasp_kpoints_setup(self):
        """Parse mast_kpoints string, which should take the format:
            number, number, number designation
            examples: "3x3x3 M", "1x1x1 G". If no designation is given,
            Monkhorst-Pack is assumed.
        """
        name = self.keywords['name']
        tryname = os.path.join(name, "KPOINTS")
        if os.path.isfile(tryname):
            #KPOINTS already exists. Do not overwrite.
            my_kpoints = Kpoints.from_file(tryname)
            return my_kpoints
        if not (self.metafile.read_data('kpoints')==None):
            kpoints = self.metafile.read_data('kpoints').split()
            kmesh = (int(kpoints[0].split('x')[0]),int(kpoints[0].split('x')[1]),int(kpoints[0].split('x')[2]))
            try: desig = kpoints[1].upper()
            except IndexError: desig = 'M'
            try: kshift = (float(kpoints[2]),float(kpoints[3]),float(kpoints[4]))
            except IndexError: kshift = (0.0,0.0,0.0)
        else:
            if 'mast_kpoints' in self.keywords['program_keys'].keys():
                kpoints = self.keywords['program_keys']['mast_kpoints']
            else:
                raise MASTError(self.__class__.__name__,"k-point instructions need to be set either in ingredients keyword mast_kpoints or scaling section in structure ingredient: No k-point settings for the ingredient %s"% name)
            if len(kpoints) == 3:
                desig = "M"
            elif 'line' in str(kpoints[1]):
                desig = 'L'
            else:
                desig = kpoints[3].upper()
            if not desig == 'L':
                kmesh = (int(kpoints[0]),int(kpoints[1]),int(kpoints[2]))
                kshift = (0,0,0)
        if desig == "M":
            my_kpoints = Kpoints.monkhorst_automatic(kpts=kmesh,shift=kshift)
        elif desig == "G":
            my_kpoints = Kpoints.gamma_automatic(kpts=kmesh,shift=kshift)
        elif desig == 'L':
            my_kpoints='Line Mode\n'+'\n'.join(' '.join(kpoints).split(','))+'\n'
            fp=open(name+'/KPOINTS','w')
            fp.write(my_kpoints)
        else:
            raise MASTError(self.__class__.__name__,"kpoint designation " + desig + " not recognized.")

        dirutil.lock_directory(name)
        if not desig=='L':
            my_kpoints.write_file(name + "/KPOINTS")
        dirutil.unlock_directory(name)
        return my_kpoints

    def _vasp_potcar_setup(self, my_poscar):
        """Set up the POTCAR file."""
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
        allowedpath = os.path.join(dirutil.get_mast_install_path(),
                        'ingredients',
                        'programkeys','vasp_allowed_keywords.py')
        allowed_list = self._vasp_incar_get_allowed_keywords(allowedpath)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                keytry = key.upper()
                if not (keytry in allowed_list):
                    self.logger.warning("Ignoring program key %s for INCAR. To allow this keyword, add it to %s" % (keytry, allowedpath))
                else:
                    if type(value)==str and value.isalpha():
                        incar_dict[keytry]=value.capitalize() #First letter cap
                    else:
                        incar_dict[keytry]=value
        return incar_dict


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
            myd['ENCUT']=self._get_max_enmax_from_potcar(my_potcar)*mymult
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
            myelectrons = self.get_total_electrons(my_poscar, my_potcar)
            newelectrons=0.0
            try:
                adjustment = float(self.keywords['program_keys']['mast_charge'])
            except (ValueError, TypeError):
                raise MASTError("vasp_checker, vasp_incar_setup","Could not parse adjustment")
            #newelectrons = myelectrons + adjustment
            newelectrons = myelectrons - adjustment
            myd['NELECT']=str(newelectrons)
        if self.metafile.read_data('nbands'):
            myd['NBANDS']=self.metafile.read_data('nbands')
        my_incar = Incar(myd)
        dirutil.lock_directory(name)
        my_incar.write_file(name + "/INCAR")
        dirutil.unlock_directory(name)
        return my_incar

    def set_up_program_input(self):
        """Set up the program input files."""
        myposcar = self._vasp_poscar_setup()
        self._vasp_kpoints_setup()
        mypotcar = self._vasp_potcar_setup(myposcar)
        self._vasp_incar_setup(mypotcar, myposcar)
        return

    def forward_extra_restart_files(self, childpath):
        """Forward extra restart files: 
            For VASP, this entails a softlink to WAVECAR and 
            CHGCAR.
            Args:
                childpath <str>: path to child ingredient
        """
        raise MASTError(self.__class__.__name__,"This function is obsolete. Use softlink_charge_density_file, softlink_wavefunction_file, forward_charge_density_file, and/or forward_wavefunction_file")
        parentpath = self.keywords['name']
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

    def softlink_charge_density_file(self, childpath):
        """Softlink the parent charge density file to the child folder
            For VASP, this entails a softlink to CHGCAR
            Args:
                childpath <str>: path to child ingredient
        """
        return self.softlink_a_file(childpath, "CHGCAR")
    
    def softlink_wavefunction_file(self, childpath):
        """Softlink the parent wavefunction file to the child folder
            For VASP, this entails a softlink to WAVECAR.
            Args:
                childpath <str>: path to child ingredient
        """
        return self.softlink_a_file(childpath, "WAVECAR")

    def forward_charge_density_file(self, childpath):
        """Copy the charge density file to the child folder.
            For VASP, this is the CHGCAR.
            Args:
                childpath <str>: path to child ingredient
        """
        return self.copy_a_file(childpath, "CHGCAR","CHGCAR")
    def forward_wavefunction_file(self, childpath):
        """Copy the wavefunction file to the child folder.
            For VASP, this is the WAVECAR.
            CHGCAR.
            Args:
                childpath <str>: path to child ingredient
        """
        return self.copy_a_file(childpath, "WAVECAR","WAVECAR")
    
    def add_selective_dynamics_to_structure_file(self, sdarray):
        """Add selective dynamics to a structure file.
            In the case of VASP, the structure file is 
            POSCAR-like
            Args:
                sdarray <numpy array of bool>: Numpy array of
                    booleans for T F F etc.
        """
        name = self.keywords['name']
        pname = os.path.join(name,"POSCAR")
        oldposcar = Poscar.from_file(pname)
        phposcar = Poscar(oldposcar.structure, selective_dynamics=sdarray.tolist())
        #phposcar.selective_dynamics = sdarray
        dirutil.lock_directory(name)
        os.rename(pname, pname + "_no_sd")
        self.write_poscar_with_zero_velocities(phposcar, pname)
        dirutil.unlock_directory(name)
        return

    def get_energy(self):
        """Get the energy.
            For VASP, this is E0 energy from vasprun.xml
        """
        return Vasprun('%s/vasprun.xml' % self.keywords['name']).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

    def _get_max_enmax_from_potcar(self, mypotcar):
        """Get maximum enmax value (float) from Potcar 
            (combined list)
        """
        enmax_list=list()
        potcarct=0
        onepotcar=None
        while potcarct < len(mypotcar):
            onepotcar = mypotcar[potcarct] #A PotcarSingle object
            enmax_list.append(onepotcar.enmax)
            potcarct = potcarct + 1
        return max(enmax_list)

    def _make_one_unfrozen_atom_poscar(self, myposcar, natom):
        """Use selective dynamics to make a poscar with one unfrozen atom.
            myposcar = Poscar
            natom = the number of the atom to unfreeze
            Returns: Poscar (use write_file function on it).
        """
        raise MASTError(self.__class__.__name__, "This method is abandoned and should not be used.")
        mysd=np.zeros([sum(myposcar.natoms),3],bool)
        mysd[natom-1][0]=True #indexing starts at 0
        mysd[natom-1][1]=True
        mysd[natom-1][2]=True
        myposcar.selective_dynamics = mysd
        return myposcar

    def _make_one_unfrozen_direction_poscar(self, myposcar, natom, ndir):
        """Use selective dynamics to make a poscar with one unfrozen atom.
            myposcar = Poscar
            natom = the number of the atom to unfreeze
            ndir = the direction to freeze (0, 1, 2 for x, y, z)
            Returns: Poscar (use write_file function on it).
        """
        raise MASTError(self.__class__.__name__, "This method is abandoned and should not be used.")
        mysd=np.zeros([sum(myposcar.natoms),3],bool)
        mysd[natom-1][ndir]=True #indexing starts at 0
        myposcar.selective_dynamics = mysd
        return myposcar


    def read_my_dynamical_matrix_file(self, mydir="", fname="DYNMAT"):
        """Read a dynamical matrix file.
            For VASP this is DYNMAT.
            Args:
                mydir <str>: directory; use ingredient directory
                             if null is given
                fname <str>: file name (default 'DYNMAT')
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
        if mydir == "":
            mydir = self.keywords['name']
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

    def write_my_dynamical_matrix_file(self, dyndict, mydir="", fname="DYNMAT"):
        """Write a dynamical matrix file based on a dictionary.
            Args:
                dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
                mydir <str>: Directory in which to write;
                             use ingredient directory if null
                fname <str>: filename (default DYNMAT)
        """
        if mydir == "":
            mydir = self.keywords['name']
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


    def write_my_dynmat_without_disp_or_mass(self, dyndict, mydir="", fname="DYNMAT"):
        """Write a dynamical matrix file without the displacement indicators 1, 2, 3
            and without the masses line, and with first line having only
            the total number of displacements, for PHON.
            Args:
                dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
                mydir <str>: Directory in which to write; use 
                             ingredient directory if null
                fname <str>: filename (default DYNMAT)
        """
        if mydir == "":
            mydir = self.keywords['name']
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

    def read_my_displacement_file(self, mydir="", fname="XDATCAR"):
        """Read a displacement file. For VASP this is XDATCAR.
            Args:
                mydir <str>: Directory. Use ingredient directory
                             if null.
                fname <str>: Filename to read. Default is
                             XDATCAR
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
        if mydir == "":
            mydir = self.keywords['name']
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

    def write_my_displacement_file(self, xdatdict, mydir="", fname="XDATCAR"):
        """Write a displacement file.
            Args:
                xdatdict <dict>: Dictionary of XDATCAR (see read_my_xdatcar)
                mydir <str>: Directory in which to write.
                             Use ingredient directory if null.
                fname <str>: filename (default DYNMAT)
        """
        if mydir == "":
            mydir = self.keywords['name']
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

    def combine_dynamical_matrix_files(self, mydir=""):
        """Combine dynamical matrix files into one file.
            The subfiles should be in subfolders under
            the top directory.
            Args:
                mydir <str>: top directory for DYNMAT files.
                             Use ingredient directory if null.
            Returns:
                Creates a single DYNMAT file
        """
        if mydir == "":
            mydir = self.keywords['name']
        dynmatlist = dirutil.walkfiles(mydir, 2, 5, "*DYNMAT*") #start one level below
        if len(dynmatlist) == 0:
            raise MASTError("pmgextend combine_dynmats", "No DYNMATs found under " + mydir)
        totnumdisp=0
        largedyn=dict()
        largedyn['atoms'] = dict()
        for onedynmat in dynmatlist:
            dyndir = os.path.dirname(onedynmat)
            onedyn = self.read_my_dynamical_matrix_file(dyndir)
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
        self.write_my_dynamical_matrix_file(largedyn, mydir, "DYNMAT_combined")
        self.write_my_dynamical_matrix_file(largedyn, mydir, "DYNMAT")

    def combine_displacement_files(self, mydir=""):
        """Combine displacement files (here XDATCARs) into one file.
            Args:
                mydir <str>: top directory for DYNMAT files.
                             Use ingredient directory if null.
        """
        if mydir == "":
            mydir = self.keywords['name']
        largexdat=dict()
        xdatlist = dirutil.walkfiles(mydir, 2, 5, "*XDATCAR*") #start one level below
        if len(xdatlist) == 0:
            raise MASTError("pmgextend combine_displacements", "No XDATCARs found under " + mydir)
        kfgct=1 # skip config 1 until the end
        largexdat['configs']=dict()
        xdatlist.sort() #get them in order
        for onexdatmat in xdatlist:
            xdatdir = os.path.dirname(onexdatmat)
            onexdat = self.read_my_displacement_file(xdatdir)
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
        self.write_my_displacement_file(largexdat, mydir, "XDATCAR_combined")
        self.write_my_displacement_file(largexdat, mydir, "XDATCAR")
    def make_hessian(self, myposcar, mydir):
        """Combine DYNMATs into one hessian and solve for frequencies.
            myposcar = Poscar
            mydir = top directory for DYNMAT files
        """
        raise MASTError(self.__class__.__name__, "This method is abandoned and should not be used.")
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

    def get_total_electrons(self, myposcar, mypotcar):
        """Get the total number of considered electrons in the system."""
        atomlist = myposcar.natoms
        zvallist = self.get_valence_list(mypotcar)
        totzval = 0.0
        atomct = 0
        if not (len(zvallist) == len(atomlist)):
            raise MASTError(self.__class__.__name__,
                "Number of species and number of POTCARs do not match.")
        while atomct < len(atomlist):
            totzval = totzval + (atomlist[atomct] * zvallist[atomct])
            atomct = atomct + 1
        return totzval

    def get_valence_list(self, mypotcar):
        """Get list of number of valence electrons for species.
            For VASP, this is a list of zvals from POTCAR.
            Args:
                mypotcar <list of PotcarSingle objects>
            Returns:
                zval_list <list of int>: List of number of
                                        valence electrons
        """
        zval_list=list()
        potcarct=0
        onepotcar=None
        while potcarct < len(mypotcar):
            onepotcar = mypotcar[potcarct] #A PotcarSingle object
            zval_list.append(onepotcar.zval)
            potcarct = potcarct + 1
        return zval_list
    def get_energy_from_energy_file(self):
        """Get the energy from the energy file.
            For VASP, this is the last E0 energy from OSZICAR,
            and this function should be used if vasprun.xml
            is corrupted or not available.
            Args:
                mydir <str>: Directory in which to look.
            Returns:
                <float>: last E0 energy from OSZICAR
        """
        fullpath=os.path.join(self.keywords['name'], "OSZICAR")
        if not os.path.isfile(fullpath):
            raise MASTError(self.__class__.__name__, "No OSZICAR file at %s" % self.keywords['name'])
        myosz = MASTFile(fullpath)
        mye0 = myosz.get_segment_from_last_line_match("E0", "E0=","d E =")
        mye0float=""
        try:
            mye0float=float(mye0)
        except TypeError:
            self.logger.error("Failed to log energy %s" % str(mye0))
        return mye0float
    def is_started(self):
        """See if the ingredient has been started on
            the queue.
        """
        if os.path.isfile(os.path.join(self.keywords['name'],'OUTCAR')):
            return True
        else:
            return False

    def write_final_structure_file(self, mystruc):
        """Write the final structure to a file.
            For VASP, this is CONTCAR.
        """
        mycontcar = Poscar(mystruc)
        cname = os.path.join(self.keywords['name'],'CONTCAR')
        self.write_poscar_with_zero_velocities(mycontcar, cname)
        return
    def has_starting_structure_file(self):
        """Evaluate whether the ingredient has a starting
            structure file. For VASP, this is a POSCAR.
        """
        return os.path.isfile(os.path.join(self.keywords['name'], 'POSCAR'))

    def has_ending_structure_file(self):
        """Evaluate whether the ingredient has a starting
            structure file. For VASP, this is a CONTCAR.
        """
        return os.path.isfile(os.path.join(self.keywords['name'], 'CONTCAR'))
    def scale_mesh(self, newstr, density):
        """Scale the kpoint mesh.
            Args:
                newstr <structure>: new structure
                density <float>: kpoint density
            Returns:
                newkmesh <pymatgen Kpoints>: new kpoint mesh
                sets program_keys 'mast_kpoints' to new mesh
        """
        newkmesh = Kpoints.automatic_density(newstr, density)
        klist = list()
        klist = newkmesh.as_dict()['kpoints'][0]
        klist.append(newkmesh.as_dict()['generation_style'][0])
        self.keywords['program_keys']['mast_kpoints'] = klist
        return newkmesh

    def get_final_pressure(self):
        """Get the final pressure.
            For VASP, this is the last pressure line from
            the OUTCAR.
            Args:
                mydir <str>: Directory in which to look.
            Returns:
                <float>: last pressure from OUTCAR, in kB
        """
        fullpath=os.path.join(self.keywords['name'], "OUTCAR")
        if not os.path.isfile(fullpath):
            raise MASTError(self.__class__.__name__, "No OUTCAR file at %s" % self.keywords['name'])
        myoutcar = MASTFile(fullpath)
        mypress = myoutcar.get_segment_from_last_line_match("pressure", "external pressure =","kB  Pullay stress =")
        mypressfloat=""
        try:
            mypressfloat=float(mypress)
        except TypeError:
            self.logger.error("Failed to log pressure %s" % str(mypress))
        return mypressfloat
