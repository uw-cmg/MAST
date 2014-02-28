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
from ase.io import read, write

def MAST_structopt(inputfile):
    """Recipe to run the optimizer through MAST
    Args:
        inputfile = Str of filename with Optimizer parameters
    """
    Opti = Optimizer(inputfile)
    Opti.algorithm_initialize()
    #Begin main algorithm loop
    while not Opti.convergence:
        if Opti.generation==0:
            pop = []
        else:
            #Clean out previous folders from MAST run
            CustomChopIngredient.clear_ga_vasp_ingredient()
        offspring = Opti.generation_set(pop)
        # Identify the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if ind.energy==0]
        #Evaluate the individuals with invalid fitness
        Opti.output.write('\n--Evaluate Structures--\n')
        #Write structures to POSCAR files for evaluation in MAST
        fobj = open('filelist.txt','w')
        for i in range(len(invalid_ind)):
            write("POSCAR_%02d"%i)
            fobj.write(os.path.join(os.getcwd(),"POSCAR_%02d"%i)+'\n')
        fobj.close()
        #Submit list of POSCARs to Tam's script to make subfolders
        make_subfolders_from_structure_list()
        #Ideally this next function will run and return a list of files that include the final
        #structure output and the energy,pressure and other output details
        Totaloutputs = CustomChopIngredient.evaluate_ga_vasp_and_update()
        #Extract information from output
        for structurefile, VASPoutfile in Totaloutputs:
            ind = read(structurefile)
            index = int(structurefile.split('_')[1])
            if Opti.structure == 'Defect':
                outt = find_defects([Opti,ind,Opti.sf,False])
                indi = outt[0].copy()
                bulki = outt[1].copy()
                invalid_ind[index][0] = indi
                invalid_ind[index].bulki = bulki
            else:
                invalid_ind[index][0] = ind
            outfile = open(VASPoutfile,'r')
            invalid_ind[index].energy = outfile.readline.split('=')
            invalid_ind[index].pressure = outfile.readline.split('=')
            # Write any other information from the structure output file to the primary 
            # optimizer output file.  This might be nothing or it could be a way to pass errors
            # through to the user without needing to preserve VASP outputs
            while True:
                try:
                    Opti.output.write(outfile.readline()+'\n')
                except:
                    break
            Opti.output.flush()
            outfile.close()
        for ind in invalid_ind:
            outs = tools.fitnesseval([Opti,ind])
            Opti.output.write(outs[1])
            invalid_ind[i]=outs[0]
        pop.extend(invalid_ind)
        pop = Opti.generation_eval(pop)
        Opti.population = pop
    
    end_signal = Opti.algorithm_stats(pop)
    return end_signal    

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
        inputfile = os.path.join(self.keywords['name'],"input.txt")
        #pathtooutput = os.path.join(self.keywords['name'],self.keywords['program_keys']['filename']) #usually filename is "GAoutput"
        pathtooutput = os.path.join(self.keywords['name'],"GAoutput")
        os.chdir(self.keywords['name'])
        Opti = Optimizer(inputfile)
        Opti.algorithm_initialize()
        #Begin main algorithm loop
        self.logger.info("opti initialized")
        genfile = MASTFile("GENERATION")
        if not os.path.isfile("GENERATION"):
            Opti.generation = 0
            pop=[]
        else:
            Opti.generation = int(genfile.data[0])
            #If VASP subfolders are not complete, return out
            #since cannot do any more checking.
            if not self.vasp_subfolders_complete():
                return 
            else: 
                #VASP subfolders are complete. Now need to
                #evaluate them for fitness and restart the
                #Genetic Algorithm loop if necessary.
                Totaloutputs = self.extract_vasp_output()
                #Extract information from output
                for (structurefile, energy, pressure) in Totaloutputs:
                    ind = read(structurefile)
                    index = int(structurefile.split('_')[1])
                    if Opti.structure == 'Defect':
                        outt = find_defects([Opti,ind,Opti.sf,False])
                        indi = outt[0].copy()
                        bulki = outt[1].copy()
                        invalid_ind[index][0] = indi
                        invalid_ind[index].bulki = bulki
                    else:
                        invalid_ind[index][0] = ind
                    invalid_ind[index].energy = energy
                    invalid_ind[index].pressure = pressure
                    # Write any other information from the structure output file to the primary 
                    # optimizer output file.  This might be nothing or it could be a way to pass errors
                    # through to the user without needing to preserve VASP outputs
                    while True:
                        try:
                            Opti.output.write(outfile.readline()+'\n')
                        except:
                            break
                    Opti.output.flush()
                    outfile.close()
                for ind in invalid_ind:
                    outs = tools.fitnesseval([Opti,ind])
                    Opti.output.write(outs[1])
                    invalid_ind[i]=outs[0]
                pop.extend(invalid_ind)
                pop = Opti.generation_eval(pop)
                Opti.population = pop
                if Opti.convergence:
                    convokay = open("CONVERGED","wb")
                    convokay.write("Opti.convergence is True at %s\n" % time.asctime())
                    convokay.close()
                    end_signal = Opti.algorithm_stats(pop)
                    return end_signal
                else:
                    self.clear_vasp_folders()
            
        offspring = Opti.generation_set(pop)
        self.logger.info("generation set")
        self.logger.info(offspring)
        # Identify the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if ind.energy==0]
        #Evaluate the individuals with invalid fitness
        Opti.output.write('\n--Evaluate Structures--\n')
        #Write structures to POSCAR files for evaluation in MAST
        fobj = open('filelist.txt','w')
        for i in range(len(invalid_ind)):
            indname = "indiv%02d.xyz" % i
            Opti.output.write(indname)
            fobj.write(os.path.join(pathtooutput,indname))
        fobj.close()
        #Submit list of POSCARs to Tam's script to make subfolders
        self.make_subfolders_from_structure_list('filelist.txt')
        #Ideally this next function will run and return a list of files that include the final
        #structure output and the energy,pressure and other output details
        self.set_up_vasp_folders()  #Start running.
        return

    def extract_vasp_output(self, childname=""):
        """Extract output from the VASP runs.
        """
        outputlist=list()
        for subfolder in self.get_vasp_subfolder_list():
            mychecker = VaspChecker(subfolder, self.keywords['program_keys'], self.keywords['structure'])
            #mystructure = mychecker.get_final_structure_from_directory()
            #structurelist.append(mystructure)
            mystrfile = os.path.join(subfolder,"CONTCAR")
            myenergy = mychecker.get_energy_from_energy_file()
            mypressure = mychecker.get_pressure()
            outputlist.append(mystrfile, myenergy, mypressure)
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

            self.change_my_status("W")

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
        for subfolder in self.get_vasp_subfolder_list():
            shutil.rmtree(subfolder)
        cleared_ing = BaseIngredient(name=fullpath)
        cleared_ing.change_my_status("W")
        return


    def make_subfolders_from_structure_list(self, fname="", startat=0, childdir=""):
        """This method makes POSCAR-containing subfolders from a 
            list of structures,
            in childdir if a child directory is given, or in 
            mydir otherwise.
            File fname must reside in directory mydir.
            Args:
                #mydir <str>: Ingredient directory (from self)
                fname <str>: File containing list of structures
                startat <str, will be converted to int>: 
                    0 (default) - start subfolders at 00
                    1 - start subfolders at 01
                childdir <str>: Child directory (optional)
            File list should be of the format:
                POSCAR_01
                POSCAR_02
                etc.
        """
        mydir = self.keywords['name']
        logger = logging.getLogger(mydir)
        logger = loggerutils.add_handler_for_recipe(mydir, logger)

        if childdir == "":
            childdir = mydir
        fpath = os.path.join(mydir, fname)
        if not os.path.isfile(fpath):
            logger.error("No file found at %s" % fpath)
            return "File not found. No effect."
        str_list_file = MASTFile(fpath)
        namelist = list(str_list_file.data)
        structure_list = self.get_structures(namelist)
        if not (len(structure_list) == len(namelist)):
            logger.error("Not all structures in file %s found." % fname)
            return "Not all structures found. No effect."
        strct = int(startat)
        for one_structure in structure_list:
            subname = str(strct).zfill(2)
            subpath = os.path.join(childdir, subname)
            if not os.path.isdir(subpath):
                os.mkdir(subpath)
            one_poscar = pymatgen.io.vaspio.Poscar(one_structure)
            pospath = os.path.join(subpath, "POSCAR")
            if os.path.isfile(pospath):
                logger.error("POSCAR file already exists at %s" % subpath)
            else:
                one_poscar.write_file(os.path.join(subpath, "POSCAR"))
            strct = strct + 1
        return "Wrote all POSCAR files to subdirectories in %s" % childdir
    def get_structures(self, namelist):
        """Get structures from list
            Args:
                #mydir <str>: Ingredient directory (from self)
                namelist <list: List of file names of structures
            Returns:
                return_structures <list of Structure objects>
        """
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
