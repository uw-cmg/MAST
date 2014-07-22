##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Amy Kaczmarowski
# Last updated: 2014-04-25
##############################################################
from pymatgen.core.structure import Structure
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.ingredients.checker import BaseChecker
import os
import logging
import pymatgen
import numpy as np
import time
import subprocess

try:
    import ase
    from ase import Atom, Atoms
except ImportError:
    print "NOTE: ASE is not installed. To use the LAMMPS checker, you must install ASE."
import random
from pymatgen.io.aseio import AseAtomsAdaptor
from re import compile as re_compile, IGNORECASE
import decimal as dec
import shutil

class LammpsChecker(BaseChecker):
    """LAMMPS checker functions
        Mostly structure functions right now.
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)
        self.potential_file = None

    def get_structure_from_file(self, myfilepath=""):
        """Get the structure from a specified file path.
            For LAMMPS, this is a trajectory-type file.
            Args:
                myfilepath <str>: File path for structure.
        """
        dir = os.path.dirname(myfilepath)
        if dir == '':
            dir = os.getcwd()
        filelist = os.listdir(dir)
        atomsymbolsfile = None
        for onefile in filelist:
            if 'atom_symbols' in onefile:
                atomsymbolsfile = os.path.join(dir, onefile)
        #Read the atom types from the atoms symbol file
        f = open(atomsymbolsfile, 'r')
        atomslist = list()
        for line in f.readlines():
            atomslist.append(line.strip())
        f.close()
        if "DATA" in myfilepath:
            atms = read_lammps_data(myfilepath,atomslist)
        elif "TRAJECTORY" in myfilepath:
            atms,velocities,forces = read_lammps_trajectory(myfilepath,-1,atomslist)
        else:
            raise MASTError(self.__class__.__name__,
                "Unknown file type for receiving structure from LAMMPS: %s" % myfilepath)
        return AseAtomsAdaptor.get_structure(atms)

    def get_initial_structure_from_directory(self,mydir=""):
        """Get the structure from a specified file path.
            For LAMMPS, this is a data-type file.
            Args:
                myfilepath <str>: File path for structure.
        """
        dir = os.path.dirname(myfilepath)
        filelist = os.listdir(dir)
        atomsymbolsfile = None
        for onefile in filelist:
            if 'atom_symbols' in onefile:
                atomsymbolsfile = os.path.join(dir, onefile)
        #Read the atom types from the atoms symbol file
        f = open(atomsymbolsfile)
        atomslist = list()
        for line in f.readlines():
            atomslist.append(line.strip())
        f.close()
        atms = read_lammps_data(myfilepath,atomslist)
        return AseAtomsAdaptor.get_structure(atms)
    
    def get_final_structure_from_directory(self, mydir=""):
        """Get the structure from a specified file path.
            For LAMMPS, this is a trajectory-type file.
            Args:
                myfilepath <str>: File path for structure.
        """
        if mydir == "":
            mydir = self.keywords['name']
        filelist = os.listdir(mydir)
        atomsymbolsfile = None
        for onefile in filelist:
            if 'atom_symbols' in onefile:
                atomsymbolsfile = os.path.join(mydir, onefile)
        #Read the atom types from the atoms symbol file
        f = open(atomsymbolsfile)
        atomslist = list()
        for line in f.readlines():
            atomslist.append(line.strip())
        f.close()
        atms = read_lammps_trajectory(myfilepath,atomslist)
        return AseAtomsAdaptor.get_structure(atms)
    
    def forward_final_structure_file(self, childpath, newname="DATA"):
        """Forward the final structure.
            For LAMMPS, this is the Trajectory.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'POSCAR')
        """
        trajstruct = self.get_structure_from_file("TRAJECTORY")
        trajatoms = AseAtomsAdaptor.get_atoms(trajstruct)
        write_lammps_data("FINAL",trajatoms)
        self.copy_a_file(childpath, "FINAL", newname)
        self.copy_a_file(childpath, "atom_symbols", "atom_symbols")
    
    def forward_initial_structure_file(self, childpath, newname="POSCAR"):
        """Forward the initial structure.
            For LAMMPS, this is the data file. This function is
            used after phonon calculations, where the CONTCAR
            contains the last displacement. To forward to PHON,
            the POSCAR (without displacements) should be used.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'POSCAR')
        """
        filelist = os.listdir(childpath)
        dataname = None
        atomsymbolsfile = None
        for onefile in filelist:
            if 'DATA' in onefile:
                dataname = os.path.join(childpath, onefile)
            if 'atom_symbols' in onefile:
                atomsymbolsfile = os.path.join(childpath, onefile)
        if not dataname:
            raise MASTError(self.__class__.__name__,"No trajectory file in %s" % mydir)
        #Read the atom types from the atoms symbol file
        f = open(atomsymbolsfile)
        atomslist = list()
        for line in f.readlines():
            atomslist.append(line.strip())
        f.close()
        atms = read_lammps_data(dataname, atomlist = atomslist)
        from ase.io import read, write
        newpath = os.path.join(childpath, "POSCAR")
        ase.io.write(newpath,aseatomsobj,"vasp", direct=True, sort=True, vasp5=True)
        return self.copy_a_file(childpath, "POSCAR", newname)
    
    def forward_displacement_file(self, childpath, newname="TRAJECTORY"):
        """Forward displacement information.
            For LAMMPS, this is the TRAJECTORY file.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'TRAJECTORY')
        """
        return self.copy_a_file(childpath, "TRAJECTORY", newname)

    def forward_energy_file(self, childpath, newname="LOG"):
        """Forward the energy file.
            For LAMMPS, this is the LOG file.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'LOG')
        """
        return self.copy_a_file(childpath, "LOG", newname)

    def is_frozen(self):
        """Check if single LAMMPS non-NEB calculation is frozen.
        """
        return BaseChecker.is_frozen(self, "LOG")
    
    def is_complete(self):
        """Check if single LAMMS non-NEB calculation is complete.
        """
        opath = os.path.join(self.keywords['name'],"LOG")
        if not os.path.isfile(opath):
            self.logger.info("No LOG at %s; not complete." % opath)
            return False
        reachgrep=subprocess.Popen('grep "ERROR" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            founderror=False
        else:
            founderror=True
            raise MASTError(self.__class__.__name__,"Found error in LAMMPS Execution: {0}".format(reachrpt))
        reachgrep=subprocess.Popen('grep "CALCULATION HAS FINISHED" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            reachedend=False
        else:
            reachedend=True
        if reachedend:
            self.logger.info("LOG file at %s shows calculation has finished." % opath)
            return True
        else:
            self.logger.info("LOG file at %s does not show that calculation has finished; still incomplete." % opath)
            return False
    
    def is_ready_to_run(self):
        """Check if single LAMMPS non-NEB ingredient is 
            ready to run.
        """
        dirname = self.keywords['name']
        flist = os.listdir(dirname)
        notready = 0
        if "DATA" not in flist:
            notready += 1
        if "INPUT" not in flist:
            notready += 1
        if self.potential_file:
            potfile = os.path.basename(self.potential_file)
            if potfile not in flist:
                notready += 1
        if "submit.sh" not in flist:
            notready = notready + 1
        if notready > 0:
            return False
        else:
            exec_line = self.keywords['program_keys']['mast_exec']
            sp_exec_line = exec_line.split(' ')
            new_exec_line = str()
            for segment in sp_exec_line:
                if ('LOG' not in segment.upper()) and ('INPUT' not in segment):
                    new_exec_line += segment + ' '
            log_name = os.path.join(dirname,"LOG")
            in_name = os.path.join(dirname,"INPUT")
            new_exec_line += "-log {0} <{1}".format(log_name,in_name)
            self.keywords['program_keys']['mast_exec'] = new_exec_line
            return True
    
    def _lammps_data_setup(self):
        """Set up the DATA file for a single LAMMPS run.
        """
        name = self.keywords['name']
        datapath = os.path.join(name, "DATA")
        if os.path.isfile(datapath):
            my_data_structure = self.get_structure_from_file(datapath)
            #parent should have given a structure
        else: #this is an originating run; mast should give it a structure
            my_data_structure = self.keywords['structure']
            self.logger.info("No DATA file found from a parent; base structure used for %s" % self.keywords['name'])
        if 'mast_coordinates' in self.keywords['program_keys'].keys():
            sxtend = StructureExtensions(struc_work1=my_data_structure, name=self.keywords['name'])
            coordstrucs=self.get_coordinates_only_structure_from_input()
            newstruc = sxtend.graft_coordinates_onto_structure(coordstrucs[0])
            my_data_structure=newstruc.copy()
        #dirutil.lock_directory(name)
        my_data_atoms = AseAtomsAdaptor.get_atoms(my_data_structure)
        write_lammps_data(datapath, my_data_atoms)
        #dirutil.unlock_directory(name)
        return datapath
    
    def _lammps_input_get_non_mast_keywords(self):
        """Get the non-LAMMPS keywords and make a dictionary."""
        input_dict=dict()
        allowedpath = os.path.join(dirutil.get_mast_install_path(),
                        'ingredients','programkeys','lammps_allowed_keywords.py')
        allowed_list = self._lammps_input_get_allowed_keywords(allowedpath)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                keytry = key.upper()
                if not (keytry in allowed_list):
                    self.logger.warning("Ignoring program key %s for INPUT. To allow this keyword, add it to %s" % (keytry, allowedpath))
                else:
                    if type(value)==str and value.isalpha():
                        input_dict[keytry.lower()]=value
                    else:
                        input_dict[keytry.lower()]=value
        return input_dict

    def _lammps_input_get_allowed_keywords(self, allowedpath):
        """Get allowed vasp keywords.
            Args:
                allowedpath <str>: file path for allowed lammps keywords
        """
        allowed = MASTFile(allowedpath)
        allowed_list=list()
        for entry in allowed.data:
            allowed_list.append(entry.strip())
        return allowed_list
    
    def _lammps_input_setup(self, lammps_input_file="INPUT", lammps_data_file="DATA", lammps_trajectory_file="TRAJECTORY"):
        """Function to write a lammps input file
        Input:
            lammps_input_file = file to write lammps input to
            lammps_data_file = file to load lammps data from
            lammps_trajectory_file = file to write lammps trajectory to
        """
        if isinstance(lammps_input_file, str):
            f = open(lammps_input_file, 'w')
        else:
            f = lammps_input_file
        
        # Write file variable names
        f.write('clear\n' +
                ('variable dump_file string "%s"\n' % lammps_trajectory_file) +
                ('variable data_file string "%s"\n' % lammps_data_file))
        f.write('units metal \n')
        
        #Write to specified log file
        #logname = os.path.join(os.path.dirname(lammps_data_file),"LOG")
        #f.write('log %s append\n' % logname)
        
        #Get keyword parameters
        parameters = dict()
        parameters = self._lammps_input_get_non_mast_keywords()
        if 'autocorrect' in parameters:
            parameters['autocorrect'] = bool(parameters['autocorrect'])
        else:
            parameters['autocorrect'] = True
        if parameters['autocorrect']:
            try:
                atpath = os.path.join(os.path.dirname(lammps_data_file),'atom_symbols')
                try:
                    atf = open(atpath,'r')
                except:
                    atf = open('atom_symbols', 'r')
                atomlist = list()
                for line in atf.readlines():
                    atomlist.append(line.strip())
                atf.close()
                atomlist = sorted(list(set(atomlist)))
            except:
                struct = self.keywords['structure']
                atomst = AseAtomsAdaptor.get_atoms(struct)
                atomlist = sorted(list(set(atomst.get_chemical_symbols())))
            parameters = readable_lammps_parameters(parameters,atomlist)
        
        #Copy potential file to current working directory
        if ('potential_file' in parameters) and (parameters['potential_file'] != 'None'):
            newlocation = os.path.join(os.path.dirname(lammps_data_file),os.path.basename(parameters['potential_file']))
            shutil.copyfile(parameters['potential_file'],newlocation)
            self.potential_file = newlocation
        
        atoms = AseAtomsAdaptor.get_atoms(self.keywords['structure'])
        pbc = atoms.get_pbc()
        if ('boundary' in parameters) and (parameters['boundary'] != 'None'):
            f.write('boundary %s \n' % parameters['boundary'])
        else:
            f.write('boundary %c %c %c \n' % tuple('sp'[x] for x in pbc))
        f.write('atom_modify sort 0 0.0 \n')
        if ('neighbor' in parameters) and (parameters['neighbor'] != 'None'):
            f.write('neighbor %s \n' % parameters['neighbor'])
        if ('newton' in parameters) and (parameters['newton'] != 'None'):
            f.write('newton %s \n' % parameters['newton'])
                
        f.write('\n')
        
        f.write('read_data %s\n' % lammps_data_file)
        
        # Write pair potential information
        f.write('\n### interactions \n')
        if ( ('pair_style' in parameters) and ('pair_coeff' in parameters) ):
            if parameters['pair_style'] == 'None':
                raise MASTError(self.__class__.__name__,"No pair_style provided for LAMMPS Input file")
            pair_style = parameters['pair_style']
            f.write('pair_style %s \n' % pair_style)
            if ('None' in parameters['pair_coeff']):
                raise MASTError(self.__class__.__name__,"No pair_coeff provided for LAMMPS Input file")
            if not isinstance(parameters['pair_coeff'],list):
                try:
                    parameters['pair_coeff'] = eval(parameters['pair_coeff'])
                except:
                    parameters['pair_coeff'] = [ parameters['pair_coeff'] ]
            for pair_coeff in parameters['pair_coeff']:
                f.write('pair_coeff %s \n' % pair_coeff)
            if ('mass' in parameters) and (parameters['mass'] != 'None'):
                if not isinstance(parameters['mass'],list):
                    try:
                        parameters['mass'] = eval(parameters['mass'])
                    except:
                        parameters['mass'] = [ parameters['mass'] ]
                for mass in parameters['mass']:
                    f.write('mass %s \n' % mass)
        else:# Use simple leonnard jones calculation
            f.write('pair_style lj/cut 2.5 \n' +
                    'pair_coeff * * 1 1 \n' +
                    'mass * 1.0 \n')
        #Set up min_style if appropriate
        if ('min_style' in parameters) and (parameters['min_style'] != 'None'):
            f.write('min_style %s \n' % parameters['min_style'])
        if ('min_modify' in parameters) and (parameters['min_modify'] != 'None'):
            f.write('min_modify %s \n' % parameters['min_modify'])
        
        #Write line for thermo arguments and stemps
        thermo_args = ['step', 'temp', 'press', 'cpu', 
                            'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                            'ke', 'pe', 'etotal',
                            'vol', 'lx', 'ly', 'lz', 'atoms']
        thermo_mark = ' '.join([x.capitalize() for x in thermo_args[0:3]])
        if ('thermosteps' not in parameters) or (parameters['thermosteps']=='None'):
            parameters['thermosteps'] = '1'
        
        f.write('\n### run\n' + 
                 'fix fix_nve all nve\n'+
                ('dump dump_all all custom {0} {1}'.format(parameters['thermosteps'], lammps_trajectory_file))+
                ' id type x y z vx vy vz fx fy fz'+
                '\n')
        f.write(('thermo_style custom %s\n' +
                'thermo_modify flush yes\n' +
                'thermo '+repr(parameters['thermosteps'])+
                '\n') % (' '.join(thermo_args)))

        rflag = False
        if ('minimize' in parameters) and (parameters['minimize'] != 'None'):
            f.write('minimize %s\n' % parameters['minimize'])
            rflag = True
        if ('run' in parameters) and (parameters['run'] != 'None'):
            f.write('run %s\n' % parameters['run'])
            rflag = True
        if not rflag:
            f.write('run 0\n')
        
        f.write('print "CALCULATION HAS FINISHED"\n')
        f.write('log /dev/stdout\n') # Force LAMMPS to flush log
        f.write('undump dump_all\n') # Force LAMMPS to flush trj
        f.close()
        return lammps_input_file
    
    def set_up_program_input(self):
        """Set up the program input files."""
        datapath = self._lammps_data_setup()
        in_name = os.path.join(os.path.dirname(datapath),"INPUT")
        traj_name = os.path.join(os.path.dirname(datapath),"TRAJECTORY")
        log_name = os.path.join(os.path.dirname(datapath),"LOG")
        self._lammps_input_setup(lammps_input_file = in_name, lammps_data_file = datapath, lammps_trajectory_file=traj_name)
        exec_line = self.keywords['program_keys']['mast_exec']
        sp_exec_line = exec_line.split(' ')
        new_exec_line = str()
        for segment in sp_exec_line:
            if ('LOG' not in segment.upper()) and ('INPUT' not in segment):
                new_exec_line += segment + ' '
        new_exec_line += "-log {0} <{1}".format(log_name,in_name)
        self.keywords['program_keys']['mast_exec'] = new_exec_line
        return

    def get_energy(self):
        """Get the energy.
        """
        thermo_content = read_lammps_log()
        return thermo_content[-1]['pe']
    
    def get_energy_from_energy_file(self):
        """Get the energy from the energy file.
            Args:
                mydir <str>: Directory in which to look.
            Returns:
                <float>: last energy from LOG file
        """
        fullpath=os.path.join(self.keywords['name'], "LOG")
        if not os.path.isfile(fullpath):
            raise MASTError(self.__class__.__name__, "No LOG file at %s" % self.keywords['name'])
        thermo_content = read_lammps_log(fullpath)
        if len(thermo_content) ==0:
            self.logger.error("Failed to log energy %s" % str(fullpath))
            return 1000
        else:
            return thermo_content[-1]['pe']

    def is_started(self):
        """See if the ingredient has been started on
            the queue.
        """
        if os.path.isfile(os.path.join(self.keywords['name'],'TRAJECTORY')):
            return True
        else:
            return False

    def write_final_structure_file(self, mystruc):
        """Write the final structure to a file.
            For LAMMPS, this is Trajectory.
        """
        write_lammps_trajectory(os.path.join(self.keywords['name'],'TRAJECTORY'),mystruc)

    def has_starting_structure_file(self):
        """Evaluate whether the ingredient has a starting
            structure file. For LAMMPS, this is a DATA file.
        """
        return os.path.isfile(os.path.join(self.keywords['name'], 'DATA'))

    def has_ending_structure_file(self):
        """Evaluate whether the ingredient has a starting
            structure file. For LAMMPS, this is a TRAJECTORY.
        """
        return os.path.isfile(os.path.join(self.keywords['name'], 'TRAJECTORY'))

    def get_final_pressure(self):
        """Get the current calculated pressure, assume isotropic medium.
            in Bar
        """
        fullpath=os.path.join(self.keywords['name'], "LOG")
        if not os.path.isfile(fullpath):
            raise MASTError(self.__class__.__name__, "No LOG file at %s" % self.keywords['name'])
        thermo_content = read_lammps_log(fullpath)
        tc = thermo_content[-1]
        stress = np.array([tc[i] for i in ('pxx','pyy','pzz')])
        return (-(stress[0] + stress[1] + stress[2]) / 3.0)

def write_lammps_data(data_name, atoms):
    """Function to write ase atoms object to LAMMPS file
    Input:
        data_name = String or fileobject for data to be written to
        atoms = ase atoms object to be written
    Output:
        None. Data file will be written
    """
    if isinstance(data_name, str):
        f = open(data_name, 'w')
    else:
        f = data_name
    f.write('\n\n')
    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    f.write('{0} \t atoms \n'.format(n_atoms))
    #Order species alphabetically
    ordered_symbols = sorted(list(set(symbols)))
    #Write these symbols to a file to ensure they are kept
    dir = os.path.dirname(data_name)
    afilename = os.path.join(dir, 'atom_symbols')
    a = open(afilename,'w')
    for sym in ordered_symbols:
        a.write('{0}\n'.format(sym))
    a.close()
    n_types = len(ordered_symbols)
    f.write('{0}  atom types\n'.format(n_types))
    p = prism(atoms.get_cell())
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()

    f.write('0.0 %s  xlo xhi\n' % xhi)
    f.write('0.0 %s  ylo yhi\n' % yhi)
    f.write('0.0 %s  zlo zhi\n' % zhi)

    if p.is_skewed():
        f.write('%s %s %s  xy xz yz\n' % (xy, xz, yz))
    f.write('\n\n')

    f.write('Atoms \n\n')
    for i, r in enumerate(map(p.pos_to_lammps_str,
                              atoms.get_positions())):
        s = ordered_symbols.index(symbols[i]) + 1
        f.write('%6d %3d %s %s %s\n' % ((i+1, s)+tuple(r)))
    f.close()

def read_lammps_trajectory(filename, timestep=-1, atomlist=False,
    ratomlist=False, writefile=False):
    """Function to convert LAMMPS Trajectory file to ASE atoms object
    Inputs:
        filename = string of LAMMPS Trajectory file to read
        timestep = Which timestep to read from the trajectory file.
            Default is last time step available
            If timestep is set to all then outputs become lists
        atomlist = list of element strings corresponding to number type
        ratomlist = Boolean for whether or not to return the atom list
        write_file = Boolean for whether or not to write ASE Atoms object to file
    Ouput:
        ASE atoms object containing structure from timestep specified
        if ratomlist: list containing element strings applied to structure corresponding to LAMMPS number
        vlist = list of velocities for each atom in structure at timestep
        flist = list of forcer for each atom in structure at timestep
    """
    f = open(filename,'r')
    if atomlist==False:
        atnum = []
        attype = []
    alist = []
    vlist = []
    flist = []
    timeflag = False
    numberflag = False
    boxflag = False
    atomsflag = False
    for line in f.readlines():
        sp = line.split()
        if timeflag == True:
            timeflag = False
        elif numberflag == True:
            noa = float(sp[0])
            numberflag = False
        elif boxflag == True:
            if boxcount <= 2:
                box.append(float(sp[1]))
                boxcount += 1
            if boxcount > 2:
                boxflag = False
        elif atomsflag == True:
            if na < noa:
                if atomlist:
                    sym = atomlist[[i for i,sym in enumerate(atomlist) if str(i+1)==sp[1]][0]]
                    at = Atom(symbol=sym, position=[float(sp[2]),float(sp[3]),float(sp[4])])
                    a.append(at)
                    v = [float(sp[i]) for i in range(5,8)]
                    vs.append(v)
                    f1 = [float(sp[i]) for i in range(8,11)]
                    fs.append(f1)
                    na += 1
                else:
                    if sp[1] in atnum:
                        sym = attype[[i for i,num in enumerate(atnum) if num==sp[1]][0]]
                        at=Atom(symbol=sym,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                        a.append(at)
                    else:
                        symn = random.choice(range(1,100))
                        at=Atom(symbol=symn,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                        a.append(at)
                        atnum.append(sp[1])
                        attype.append(at.symbol)
                    v = [float(sp[i]) for i in range(5,8)]
                    vs.append(v)
                    f1 = [float(sp[i]) for i in range(8,11)]
                    fs.append(f1)
                    na += 1
            if na >= noa:
                atomsflag = False
                if atomlist==False:
                    atomlist=list(np.zeros(len(attype)))
                    for i in range(len(atnum)):
                        atomlist[int(atnum[i])-1]=attype[i]
                alist.append(a)
                vlist.append(vs)
                flist.append(fs)
        elif sp[1] == 'TIMESTEP':
            timeflag = True
        elif sp[1] == 'NUMBER':
            numberflag = True
        elif sp[1] == 'BOX':
            boxflag = True
            box = []
            boxcount = 0
            pbcc = [False,False,False]
            if sp[3] == 'pp':
                pbcc[0] = True
            if sp[4] == 'pp':
                pbcc[1] = True
            if sp[5] == 'pp':
                pbcc[2] = True
        elif sp[1] == 'ATOMS':
            atomsflag = True
            na = 0
            a = Atoms(cell=box,pbc=pbcc)
            vs = []
            fs = []
    f.close()
    if timestep == 'all':
        writefile=False
    if writefile == True:
        write_xyz(filename+'.xyz',alist[timestep],'xyz')
    if timestep == 'all':
        if ratomlist:
            return alist, atomlist, vlist, flist
        else:
            return alist, vlist, flist
    else:
        if ratomlist:
            return alist[timestep], atomlist, vlist[timestep], flist[timestep]
        else:
            return alist[timestep], vlist[timestep], flist[timestep]

def write_lammps_trajectory(filename, structure):
    if isinstance(filename,str):
        f = open(filename,'w')
    else:
        f = filename
    f.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
    atoms = AseAtomsAdaptor.get_atoms(structure)
    f.write("{0}\n".format(len(atoms)))
    f.write("ITEM: BOX BOUNDS ")
    pbc = atoms.get_pbc()
    for bound in pbc:
        if bound:
            f.write("pp ")
        else:
            f.write("s ")
    f.write("\n")
    p = prism(atoms.get_cell())
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()

    f.write('0.0 {0}\n'.format(xhi))
    f.write('0.0 {0}\n'.format(yhi))
    f.write('0.0 {0}\n'.format(zhi))
    if p.is_skewed():
        f.write('%s %s %s\n' % (xy, xz, yz))
    f.write('ITEM: ATOMS id type x y z vx vy vz fx fy fz\n')
    syms = sorted(list(set(atoms.get_chemical_symbols())))
    symlist = [one for one in enumerate(syms)]
    for idx in range(len(atoms)):
        type = [ty+1 for ty,sym in symlist if sym==atoms[idx].symbol][0]
        f.write('{0} {1} {2}'.format(idx+1,type, ' '.join([repr(x) for x in atoms[0].position])))
        f.write('0 0 0 0 0 0\n')
    f.close()
    return

def read_lammps_data(filename, atomlist=False, ratomlist=False, writefile=False):
    """Function to read a Lammps data file and output an xyz and atoms object
    Input:
        Filename = string for data file to read
        atomlist = List of strings corresponding to Lammps atom type
            ['Cr','Fe','He',...]
        ratomlist = True/False - will return atomlist with structure
        writefile = True/False - will output an xyz file of the Lammps data file
    """
    f=open(filename,'r')
    for i in range(4):
        ln=f.readline()
    box = []
    for i in range(3):
        ln=f.readline().split()
        box.append(float(ln[1]))
    for i in range(4):
        ln=f.readline()
    a=Atoms(cell=box)
    if atomlist:
        atlist=[pair for pair in enumerate(atomlist)]
        for line in f.readlines():
            sp=line.split()
            for i,sym in atlist:
                if sp[1]==str(i+1):
                    at=Atom(symbol=sym,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                    a.append(at)
    else:
        atnum = []
        attype = []
        for line in f.readlines():
            sp=line.split()
            if sp[1] in atnum:
                sym = attype[[i for i,num in enumerate(atnum) if num==sp[1]][0]]
                at=Atom(symbol=sym,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                a.append(at)
            else:
                symn = random.choice(range(1,100))
                at=Atom(symbol=symn,position=[float(sp[2]),float(sp[3]),float(sp[4])])
                a.append(at)
                atnum.append(sp[1])
                attype.append(at.symbol)
        atomlist=list(np.zeros(len(attype)))
        for i in range(len(atnum)):
            atomlist[int(atnum[i])-1]=attype[i]
    f.close()
    if writefile==True:
        write_xyz(filename+'.xyz',a,'xyz')
    if ratomlist:
        return a, atomlist
    else:
        return a

def read_lammps_log(lammps_log="LOG"):
    """Method which reads a LAMMPS output log file."""

    if isinstance(lammps_log, str):
        f = open(lammps_log, 'r')
    else:
        f = lammps_log
    #thermo arguments
    thermo_args = ['step', 'temp', 'press', 'cpu', 
                        'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                        'ke', 'pe', 'etotal',
                        'vol', 'lx', 'ly', 'lz', 'atoms']
    thermo_mark = ' '.join([x.capitalize() for x in thermo_args[0:3]])
    f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
    n = len(thermo_args)
    thermo_re = re_compile(r'^\s*' + r'\s+'.join([f_re]*n) + r'\s*$', flags=IGNORECASE)
    thermo_content = []
    line = f.readline()
    while line and "CALCULATION HAS FINISHED" not in line.strip():
        if line.startswith(thermo_mark):
            m = True
            while m:
                line = f.readline()
                m = thermo_re.match(line)
                if m:
                    # create a dictionary between each of the thermo_style args
                    # and it's corresponding value
                    thermo_content.append(dict(zip(thermo_args, map(float, m.groups()))))
        else:
            line = f.readline()
    f.close()
    return thermo_content

def readable_lammps_parameters(parameters, atomlist):
    #atoms = AseAtomsAdaptor.get_atoms(structure)
    #atomlist = sorted(list(set(atoms.get_chemical_symbols())))
    masslist = [Atom(sym).mass for sym in atomlist]
    if 'potential_file' in parameters:
        if parameters['potential_file']:
            parameters['pot_file'] = os.path.basename(parameters['potential_file'])
    if parameters['pair_style']=='tersoff':
        parcoff = '* * {0}'.format(parameters['pot_file'])
        for one in atomlist:
            parcoff+=' {0}'.format(one[0])
        parameters['pair_coeff'] = [parcoff]
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2,masslist[i+1]))
        parameters['mass'] = mass
    elif parameters['pair_style']=='eam':
        pair_coeff = [ '* * {0}'.format(parameters['pot_file'])]
        parameters['pair_coeff'] = pair_coeff
    elif parameters['pair_style']=='eam/fs':
        parcoff = '* * {0}'.format(parameters['pot_file'])
        for one in atomlist:
            parcoff+=' {0}'.format(one[0])
        parameters['pair_coeff'] = [parcoff]
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2,masslist[i+1]))
    elif parameters['pair_style']=='eam/cd':
        parcoff = '* * {0}'.format(parameters['pot_file'])
        for one in atomlist:
            parcoff+=' {0}'.format(one[0])
        parameters['pair_coeff'] = [parcoff]
    elif parameters['pair_style']=='edip':
        parcoff = '* * {0}'.format(parameters['pot_file'])
        for one in atomlist:
            parcoff+=' {0}'.format(one)
        parameters['pair_coeff'] = [parcoff]
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2,masslist[i+1]))
        parameters['mass'] = mass
        parameters['newton'] = 'on'
    elif parameters['pair_style']=='bop':
        parcoff = '* * {0}'.format(parameters['pot_file'])
        for one in atomlist:
            parcoff+=' {0}'.format(one[0])
        parcoff+='\ncommunicate single cutoff {0}'.format(parameters['bopcutoff'])
        parameters['pair_coeff'] = [parcoff]
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2,masslist[i+1]))
        parameters['mass'] = mass
        parameters['newton'] = 'on'
    elif parameters['pair_style']=='buck':
        pairstyle='{0} {1}'.format(parametes['pair_style'], parameters['buckcutoff'])
        pair_coeff=parameters['buckparameters']
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2, masslist[i+1]))	
        parameters['pair_style'] = pairstyle
        parameters['pair_coeff'] = pair_coeff
        parameters['mass'] = mass
        parameters['potential_file'] = 'None'
    elif parameters['pair_style']=='other':
        mass = ['1 {0}'.format(masslist[0])]
        if len(masslist) > 1:
            for i in range(len(masslist)-1):
                mass.append('{0} {1}'.format(i+2,masslist[i+1]))
        if 'charges' in parameters:
            cs = parameters['charges']
            parameters['mass'][len(parameters['mass'])-1] +=cs[1]
            parameters['newton']+='\natom_style charge'
    else:
        parameters={'potential_file':'None'}
        print 'WARNING: No LAMMPS potential recognized. Assuming Lennard Jones Potential'
    return parameters

class prism:
    def __init__(self, cell, pbc=(True,True,True), digits=10):
        """COPYRIGHT: ATOMISTIC SIMULATION ENVIRONMENT
        PLEASE SEE https://wiki.fysik.dtu.dk/ase/ FOR SIMILAR METHODS
        AND MORE INFORMATION ON THIS CLASS STRUCTURE
        
        Create a lammps-style triclinic prism object from a cell

        The main purpose of the prism-object is to create suitable 
        string representations of prism limits and atom positions
        within the prism.
        When creating the object, the digits parameter (default set to 10)
        specify the precission to use.
        lammps is picky about stuff being within semi-open intervals,
        e.g. for atom positions (when using create_atom in the in-file), 
        x must be within [xlo, xhi).
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]
        
        alpha = np.arccos(np.dot(b, c)/(bn*cn))
        beta  = np.arccos(np.dot(a, c)/(an*cn))
        gamma = np.arccos(np.dot(a, b)/(an*bn))
        
        xhi = an
        xyp = np.cos(gamma)*bn
        yhi = np.sin(gamma)*bn
        xzp = np.cos(beta)*cn
        yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    
        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(max((xhi,yhi,zhi))))-digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0,   0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5*p
            n = (np.mod(x, p) - x)/p
            return [float(self.f2qdec(a)) for a in (vec + n*pvec)]

        Apre[1,:] = fold(Apre[1,:], Apre[0,:], 0)
        Apre[2,:] = fold(Apre[2,:], Apre[1,:], 1)
        Apre[2,:] = fold(Apre[2,:], Apre[0,:], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2qs(self, f):
        return str(self.f2qdec(f))

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_HALF_EVEN))

    def dir2car(self, v):
        "Direct to cartesian coordinates"
        return np.dot(v, self.A)

    def car2dir(self, v):
        "Cartesian to direct coordinates"
        return np.dot(v, self.Ainv)

    def fold_to_str(self,v):
        "Fold a position into the lammps cell (semi open), return a tuple of str" 
        # Two-stage fold, first into box, then into semi-open interval
        # (within the given precission).
        d = [x % (1-self.dir_prec) for x in 
             map(dec.Decimal, map(repr, np.mod(self.car2dir(v) + self.eps, 1.0)))]
        return tuple([self.f2qs(x) for x in 
                      self.dir2car(map(float, d))])
        
    def get_lammps_prism(self):
        A = self.A
        return (A[0,0], A[1,1], A[2,2], A[1,0], A[2,0], A[2,1])

    def get_lammps_prism_str(self):
        "Return a tuple of strings"
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])

    def pos_to_lammps_str(self, position):
        "Rotate an ase-cell postion to the lammps cell orientation, return tuple of strs"
        return tuple([self.f2s(x) for x in np.dot(position, self.R)])

    def pos_to_lammps_fold_str(self, position):
        "Rotate and fold an ase-cell postion into the lammps cell, return tuple of strs"
        return self.fold_to_str(np.dot(position, self.R))

    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)
        
