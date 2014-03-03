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
from structopt.tools import find_defects, fitness_switch
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
        gadir = Opti.__dict__['filename'] #usually GAoutput
        pathtooutput = os.path.join(ingpath,gadir)
        perfect_structure = ase.io.read(os.path.join(ingpath, Opti.__dict__['SolidFile']))
        rcutoff=5.0
        os.chdir(self.keywords['name'])
        Opti.algorithm_initialize()
        #Begin main algorithm loop
        self.logger.info("opti initialized")
        pop=[]
        if not os.path.isfile("GENERATION"):
            self.logger.info("Generation 0")
            genfile = MASTFile()
            Opti.generation = 0
        else:
            genfile = MASTFile("GENERATION")
            Opti.generation = int(genfile.data[0])
            self.logger.info("Generation %i" % Opti.generation)
            #If VASP subfolders are not complete, return out
            #since cannot do any more checking.
            if not self.vasp_subfolders_complete():
                self.logger.info("Subfolders not complete. No further checks.")
                #Opti.output.close()
                return 
            else: 
                #VASP subfolders are complete. Now need to
                #evaluate them for fitness and restart the
                #Genetic Algorithm loop if necessary.
                self.logger.info("Extract VASP outputs.")
                Totaloutputs = self.extract_vasp_output()
                #Append the new structures onto the indiv xyz's.
                for idx in range(0, len(Totaloutput)):
                    (structurefile, energy, pressure) = Totaloutput[idx]
                    indivname = "indiv%s.xyz" % idx.zfill(2)
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
                Opti.restart = True
                offspring = Opti.generation_set(pop)
                self.logger.info("Generation set")
                self.logger.info(offspring)
                #Add energy and pressure info onto offspring.
                for idx in range(0, len(Totaloutput)):
                    (structurefile, energy, pressure) = Totaloutput[idx]
                    offspring[idx].energy = energy
                    offspring[idx].pressure = pressure
                #Set fitnesses of offspring
                for ind in offspring:
                    outs = fitness_switch([Opti,ind])
                    Opti.output.write(outs[1])
                    ind=outs[0]
                #Build up the offspring, with fitness information
                pop.extend(offspring)
                #Evaluate generation (set rank of fitnesses)
                pop = Opti.generation_eval(pop)
                Opti.population = pop
                #Get previous convergence data
                from MAST.utility import fileutil
                outputpath = "%s/%s.txt" % (pathtooutput, gadir)
                endmin = fileutil.grepme(outputpath, "Min", 20)
                endavg = fileutil.grepme(outputpath, "Avg", 20)
                endgenrep = fileutil.grepme(outputpath, "Genrep", 20)
                if Opti.convergence_scheme == "Gen-Rep-Min":
                    Opti.minfit = endmin
                elif Opti.convergence_scheme == "Gen-Rep-Max":
                    Opti.minfit = endavg
                Opti.genrep = endgenrep
                #Check convergence
                check_results = Opti.check_pop()

                if Opti.convergence:
                    self.logger.info("Converged!")
                    convokay = open("CONVERGED","wb")
                    convokay.write("Opti.convergence is True at %s\n" % time.asctime())
                    convokay.close()
                    end_signal = Opti.algorithm_stats(pop)
                    #Opti.output.close()
                    return end_signal
                else:
                    self.logger.info("Clearing folders for next setup.")
                    self.clear_vasp_folders()
            
        #Mutating and evolving
        offspring = Opti.generation_set(pop)
        self.logger.info("generation set")
        self.logger.info(offspring)
        # Identify the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if ind.energy==0]
        #Write structures to POSCAR files for evaluation in MAST
        for i in range(len(invalid_ind)):
            indname = "%s/POSCAR_%02d" % (gadir, i)
            ase.io.write(indname,invalid_ind[i][0],"vasp")
            Opti.output.write(indname)
        #Submit list of POSCARs to Tam's script to make subfolders
        self.make_subfolders_from_structures(pathtooutput)
        #Ideally this next function will run and return a list of files that include the final
        #structure output and the energy,pressure and other output details
        self.set_up_vasp_folders()  #Start running.
        generation = Opti.generation + 1
        genfile.data="%i\n" % generation
        genfile.to_file(os.path.join(ingpath, "GENERATION"))
        cleared_ing = ChopIngredient(name=ingpath)
        cleared_ing.change_my_status("S")
        #Opti.output.close()
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
        if allcomplete == len(subfolders):
            return True
        else:
            return False
