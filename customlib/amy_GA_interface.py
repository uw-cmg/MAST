import os
import shutil
import sys
import time
import logging
from MAST.utility import loggerutils
from MAST.utility import MASTFile
import numpy as np
import pymatgen
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.chopingredient import ChopIngredient
from pymatgen.core.structure import Structure
from structopt import Optimizer
from structopt.generate import Individual
from structopt.tools import find_defects, fitness_switch
from structopt import tools
import structopt

#Need the following for "eval" in fitness_switch
from structopt.fitness import *
import ase

class OptiIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        """Please not modify this init method."""
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

    def evaluate(self):
        ingpath = self.keywords['name']
        inputfile = os.path.join(ingpath,"input.txt")
        Opti = Optimizer(inputfile)
        Opti.restart_ints = 0
        Opti.cxattempts = 0
        Opti.mutattempts = []
        #self.logger.info("OPTI VARS:")
        #optivars = vars(Opti)
        #for item in optivars.items():
        #    self.logger.info("%20s: %20s" % item)
        os.chdir(self.keywords['name'])
        Opti.algorithm_initialize()
        #Begin main algorithm loop
        self.logger.info("Opti initialized")
        if not os.path.isfile("GENERATION"):
            self.logger.info("Generation 0")
            genfile = MASTFile()
            Opti.generation = 0
            Opti.population=[]
        else:
            genfile = MASTFile("GENERATION")
            Opti.generation = int(genfile.data[0]) 
            self.logger.info("Generation %i" % Opti.generation)
            #If VASP subfolders are not complete, return out
            #since cannot do any more checking.
            if not self.vasp_subfolders_complete():
                self.logger.info("Subfolders not complete. No further checks.")
                #Opti.output.close()
                self.flip_status_back()
                return 
            else: 
                self.logger.info("Evaluating completed structures.")
                Opti = self.evaluate_structures_for_fitness(Opti)
                if Opti.convergence:
                    return
        genfile.data="%i\n" % Opti.generation
        genfile.to_file(os.path.join(ingpath, "GENERATION"))
        self.logger.info("Mutating and evolving structures for generation %i" % Opti.generation)
        self.mutate_evolve_write_new_structures(Opti)
        self.logger.info("Changing status back.")
        self.flip_status_back()
        #Opti.output.close()
        return

    def evaluate_structures_for_fitness(self, MyOpti):
        """Evaluate completed VASP runs' structures for fitness.
            Args:
                MyOpti <GA Optimizer>
            Returns:
                MyOpti <GA Optimizer>: optimizer with any 
                        settings that had been set
        """
        #VASP subfolders are complete. Now need to
        #evaluate them for fitness and restart the
        #Genetic Algorithm loop if necessary.
        self.logger.info("Extracting VASP outputs.")
        ingpath = self.keywords['name']
        pathtooutput = os.path.join(ingpath,MyOpti.filename)
        Totaloutputs = self.extract_vasp_output()
        #Append the new structures onto the indiv xyz's.
        for idx in range(0, len(Totaloutputs)):
            (structurefile, energy, pressure) = Totaloutputs[idx]
            indivname = "indiv%s.xyz" % str(idx).zfill(2)
            pathtovasp = os.path.dirname(structurefile)
            atomsobj = ase.io.read(structurefile)
            ase.io.write(os.path.join(pathtovasp, indivname), atomsobj, "xyz")
            newxyz = open(os.path.join(pathtovasp,indivname),'rb')
            newlines = newxyz.readlines()
            newxyz.close()
            oldxyz = open(os.path.join(pathtooutput,indivname),'ab') #Append
            oldxyz.writelines(newlines)
            oldxyz.close()
        #Set a new generation from the structures.
        MyOpti.restart = True
        pop=[]
        offspring = MyOpti.generation_set(pop)
        self.logger.info("Generation has been set.")
        #self.logger.info(offspring)
        #Add energy and pressure info onto offspring.
        for idx in range(0, len(Totaloutputs)):
            (structurefile, energy, pressure) = Totaloutputs[idx]
            offspring[idx].energy = energy
            offspring[idx].pressure = pressure
        #Set fitnesses of offspring
        self.logger.info(os.environ)
        for ind in offspring:
            outs = fitness_switch([MyOpti,ind])
            MyOpti.output.write(outs[1])
            ind=outs[0]
        #Build up the offspring, with fitness information
        pop.extend(offspring)
        #Read in BESTS as Individuals with fitnesses
        if MyOpti.BestIndsList and (MyOpti.generation > 1):
            try:
                MyOpti.BESTS = structopt.io.read_xyz(os.path.join(pathtooutput, "Bests-%s.xyz" % MyOpti.filename),n='All',data=False)
                MyOpti.BESTS = [Individual(ind) for ind in MyOpti.BESTS]
                fitlist = open("Bests-energies%s.txt" % MyOpti.filename,'rb')
                fitlines = fitlist.readlines()
                for fdx in range(0, len(fitlines)):
                    MyOpti.BESTS[fdx].fitness = float(fitlines[fdx])
                    MyOpti.BESTS[fdx].bulki = ase.Atoms()
            except IOError:
                self.logger.error("No BESTS file found.")
                MyOpti.BESTS=[]
        if MyOpti.generation < 2:
            MyOpti.BESTS=[]
        #Get previous convergence data
        from MAST.utility import fileutil
        MyOpti.output.flush()
        outputpath = "%s/%s.txt" % (ingpath, MyOpti.filename)
        hasstats = fileutil.grepme(outputpath, "Stats")
        if not hasstats == []:
            grepmin = fileutil.grepme(outputpath, "  Min ")
            if not grepmin == []:
                endmin = float(grepmin[-1].split()[1])
            else:
                endmin = None

            grepavg = fileutil.grepme(outputpath, "  Avg ")
            if not grepavg == []:
                endavg = float(grepavg[-1].split()[1])
            else:
                endavg = None

            grepgenrep = fileutil.grepme(outputpath, "  Genrep ")
            if not grepgenrep == []:
                endgenrep = float(grepgenrep[-1].split()[1])
            else:
                endgenrep = None
            if MyOpti.convergence_scheme == "Gen_Rep_Min":
                MyOpti.minfit = endmin
            elif MyOpti.convergence_scheme == "Gen_Rep_Max":
                MyOpti.minfit = endavg
            MyOpti.genrep = endgenrep
        #Evaluate generation (set rank of fitnesses)
        pop = MyOpti.generation_eval(pop)
        MyOpti.population = pop
        #Check convergence
        #check_results = Opti.check_pop(pop)
        MyOpti.generation = MyOpti.generation + 1
        self.logger.info("Generation incremented to %i" % MyOpti.generation)
        if MyOpti.convergence:
            self.logger.info("Converged!")
            convokay = open("CONVERGED","wb")
            convokay.write("Opti.convergence is True at %s\n" % time.asctime())
            convokay.close()
            end_signal = MyOpti.algorithm_stats(pop)
            #Opti.output.close()
            return end_signal
        else:
            self.logger.info("Clearing folders for next setup.")
            self.clear_vasp_folders()
        return MyOpti

    def mutate_evolve_write_new_structures(self, MyOpti):
        """Mutate, evolve, and write new structures for VASP.
            Args:
                MyOpti <GA Optimizer>
            Returns:
                MyOpti <GA Optimizer>: optimizer with any 
                        settings that had been set
        """
        #Mutating and evolving
        ingpath = self.keywords['name']
        pathtooutput = os.path.join(ingpath,MyOpti.filename)
        offspring = MyOpti.generation_set(MyOpti.population)
        self.logger.info("generation set; pop len %i" % len(MyOpti.population))
        # Identify the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if ind.energy==0]
        #Write structures to POSCAR files for evaluation in MAST
        for i in range(len(invalid_ind)):
            indname = "%s/POSCAR_%02d_ase" % (MyOpti.filename, i)
            if MyOpti.structure == 'Defect':
                indatoms = invalid_ind[i][0]
                indatoms.extend(invalid_ind[i].bulki)
                ase.io.write(indname, indatoms, "vasp", vasp5=True)
            else:
                ase.io.write(indname,invalid_ind[i][0],"vasp", vasp5=True)
            mypos = pymatgen.io.vaspio.Poscar.from_file(indname)
            mypos.write_file("%s/POSCAR_%02d" % (MyOpti.filename, i))
            MyOpti.output.write(indname)
        #Submit list of POSCARs to Tam's script to make subfolders
        self.make_subfolders_from_structures(pathtooutput)
        #Ideally this next function will run and return a list of files that include the final
        #structure output and the energy,pressure and other output details
        self.set_up_vasp_folders()  #Start running.
        return MyOpti

    def flip_status_back(self):
        """Flip status back to waiting."""
        ingpath = self.keywords['name']
        cleared_ing = ChopIngredient(name=ingpath, program='None',program_keys = self.keywords['program_keys'], structure=self.keywords['structure'])
        cleared_ing.change_my_status("W")
        return

    def extract_vasp_output(self, childname=""):
        """Extract output from the VASP runs.
        """
        outputlist=list()
        for subfolder in self.get_vasp_subfolder_list():
            mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            #mystructure = mychecker.get_final_structure_from_directory()
            #structurelist.append(mystructure)
            mystrfile = os.path.join(subfolder,"CONTCAR")
            myenergy = mychecker.get_energy_from_energy_file()
            mypressure = mychecker.get_final_pressure()
            outputlist.append((mystrfile, myenergy, mypressure))
        return outputlist

    def set_up_vasp_folders(self, childname=""):
        """Set up the Genetic Algorithm VASP ingredient subfolders
        """
        for subfolder in self.get_vasp_subfolder_list():
            mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            mychecker.set_up_program_input()
            mychoping = ChopIngredient(name=subfolder, program='vasp',program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
            mychoping.write_submit_script()
            mychoping.run_singlerun()

    def clear_vasp_folders(self):
        """Clear the Genetic Algorithm VASP ingredient.
            Assumes that the VASP ingredient has already
            been evaluated and the result passed
            to the child ingredient.
            Args:
                ingred_name <str>: Ingredient name to clear.
                    Match this carefully with a name from
                    the recipe.
        """
        fullpath = self.keywords['name']
        self.logger.info("Removing directories and files from ingredient specified at %s" % fullpath)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        os.mkdir(timestamp)
        for subfolder in self.get_vasp_subfolder_list():
            os.rename(subfolder, os.path.join(fullpath, timestamp, os.path.basename(subfolder)))
            #shutil.rmtree(subfolder)
        return


    def make_subfolders_from_structures(self, finddir="", childdir=""):
        """This method makes POSCAR-containing subfolders from POSCAR
            files found in a directory
            in childdir if a child directory is given, or in 
            mydir otherwise.
            Args:
                #mydir <str>: Ingredient directory (from self)
                finddir <str>: Directory containing list of structures
                childdir <str>: Child directory (optional)
            Files should be named:
                POSCAR_00
                POSCAR_01
                POSCAR_02
                etc.
        """
        mydir = self.keywords['name']
        logger = logging.getLogger(mydir)
        logger = loggerutils.add_handler_for_recipe(mydir, logger)
        if childdir == "":
            childdir = mydir
        strfiles = os.listdir(finddir)
        for strfile in strfiles:
            if not strfile[0:6] == "POSCAR":
                continue
            if "ase" in strfile:
                continue
            str_path = os.path.join(finddir, strfile)
            subname = strfile.split("_")[1]
            subpath = os.path.join(childdir, subname)
            if not os.path.isdir(subpath):
                os.mkdir(subpath)
            shutil.copy(str_path, "%s/POSCAR" % subpath)
        return "Wrote all POSCAR files to subdirectories in %s" % childdir
    def get_structures(self, namelist):
        """Get structures from list
            Args:
                #mydir <str>: Ingredient directory (from self)
                namelist <list: List of file names of structures
            Returns:
                return_structures <list of Structure objects>
        """
        raise NotImplementedError()
        mydir = self.keywords['name']
        logger = logging.getLogger(mydir)
        logger = loggerutils.add_handler_for_recipe(mydir, logger)
        return_structures = list()
        
        for strline in namelist:
            tryline = strline.strip()
            if os.path.isfile(tryline):
                str_path = tryline
            else:
                trypath = os.path.join(mydir, tryline)
                if os.path.isfile(trypath):
                    str_path = trypath
                else:
                    logger.error("No file found at %s" % tryline)
                    continue
            if 'POSCAR' in str_path:
                one_structure = pymatgen.io.smartio.read_structure(str_path)
            elif 'xyz' in str_path:
                myxyz = pymatgen.io.xyzio.XYZ.from_file(str_path)
                allarr = np.array(myxyz.molecule.cart_coords)
                #Get box size:
                maxa=max(allarr[:,0])
                maxb=max(allarr[:,1])
                maxc=max(allarr[:,2])
                boxtol=0.01
                one_structure = myxyz.molecule.get_boxed_structure(maxa+boxtol,maxb+boxtol,maxc+boxtol)
            return_structures.append(one_structure)
        return return_structures
    def get_vasp_subfolder_list(self):
        """Get VASP subfolder list (numeric digits)
        """
        dircontents = os.listdir(self.keywords['name'])
        subfolders = list()
        for diritem in dircontents:
            fulldir = os.path.join(self.keywords['name'],diritem)
            if os.path.isdir(fulldir) and diritem.isdigit():
                subfolders.append(fulldir)
        subfolders.sort()
        return subfolders

    def vasp_subfolders_complete(self, childname=""):
        """Set up the Genetic Algorithm VASP ingredient subfolders
        """
        subfolders = self.get_vasp_subfolder_list()
        allcomplete=0
        for subfolder in self.get_vasp_subfolder_list():
            mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            if mychecker.is_complete():
                allcomplete = allcomplete + 1
            else:
                mychoping = ChopIngredient(name=subfolder, program='vasp',program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
                mychoping.run_singlerun()

        if allcomplete == len(subfolders) and (allcomplete > 0):
            return True
        else:
            return False
