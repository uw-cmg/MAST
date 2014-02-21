import os
import shutil
import sys
import time
import logging
from MAST.utility import loggerutils
from MAST.ingredients.baseingredient import BaseIngredient
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
        os.chdir(self.keywords['name'])
        Opti = Optimizer(inputfile)
        Opti.algorithm_initialize()
        #Begin main algorithm loop
        self.logger.info("opti initialized")
        while not Opti.convergence:
            if Opti.generation==0:
                pop = []
            else:
                #Clean out previous folders from MAST run
                self.clear_ga_vasp_ingredient()
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
                write("POSCAR_%02d"%i)
                fobj.write(os.path.join(os.getcwd(),"POSCAR_%02d"%i)+'\n')
            fobj.close()
            #Submit list of POSCARs to Tam's script to make subfolders
            self.make_subfolders_from_structure_list('filelist.txt')
            #Ideally this next function will run and return a list of files that include the final
            #structure output and the energy,pressure and other output details
            self.set_up_VASP_folders()
            Totaloutputs = self.evaluate_ga_vasp_and_update()
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

    def my_custom_chop_method(self, optional_string=""):
        self.logger.info("Here is a custom function within a custom class.")
        self.logger.info("This function can access the ingredient's program keys.")
        self.logger.info("%s" % self.keywords['program_keys'])
        return "Hello! %s" % optional_string

    def evaluate_ga_vasp_and_update(self, childname=""):
        """Evaluate the Genetic Algorithm VASP ingredient.
        """
        raise NotImplementedError
        childpath = os.path.join(os.path.dirname(self.keywords['name']), childname)
        from mastlib.amy_ga_code import fitness_evaluation
        from MAST.ingredients.checker import VaspChecker
        from MAST.utility import MASTFile
        dircontents = os.listdir(self.keywords['name'])
        subfolders = list()
        for diritem in dircontents:
            fulldir = os.path.join(self.keywords['name'],diritem)
            if os.path.isdir(fulldir) and diritem.isdigit():
                subfolders.append(fulldir)
        
        energylist = list()
        structurelist = list()
        for subfolder in subfolders:
            mychecker = VaspChecker(subfolder, self.keywords['program_keys'], self.keywords['structure'])
            mystructure = mychecker.get_final_structure_from_directory()
            structurelist.append(mystructure)
            myenergy = mychecker.get_energy_from_energy_file()
            energylist.append(myenergy)

        [fitoutput, fitstructure] = fitness_evaluation.evaluate(structurelist, energylist)
        #If output is a structure or xyz file, could just write it directly.
        fitfile = MASTFile()
        fitfile.data = fitoutput
        import time
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        outputname = "my_output_%s" % timestamp
        outputstrname = "my_structure_%s" % timestamp
        fitfile.to_file(os.path.join(childpath, outputname)) 
        fitstructure.write_file(os.path.join(childpath, outputstrname))
        return " %s and %s written in %s" % (outputname, outputstrname, childpath)

    def set_up_VASP_folders(self, childname=""):
        """Set up the Genetic Algorithm VASP ingredient subfolders
        """
        from MAST.ingredients.checker import VaspChecker
        from MAST.ingredients.chopingredient import ChopIngredient
        from MAST.utility import MASTFile
        dircontents = os.listdir(self.keywords['name'])
        subfolders = list()
        for diritem in dircontents:
            fulldir = os.path.join(self.keywords['name'],diritem)
            if os.path.isdir(fulldir) and diritem.isdigit():
                subfolders.append(fulldir)
        for subfolder in subfolders:
            mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            mychecker.set_up_program_input()
            mychoping = ChopIngredient(name=subfolder, program='vasp',program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
            mychoping.write_submit_script()

    def evaluate_ga_all_fitness_outputs(self):
        """Evaluate all fitness outputs"""
        outputname = "my_output"
        dircontents = os.listdir(self.keywords['name'])
        filelist = list()
        for myfile in dircontents:
            if outputname in myfile:
                filelist.append(os.path.join(self.keywords['name'], myfile))
        filelist.sort()
        last_file = filelist[-1]
        second_to_last_file = filelist[-2]
        from amy_ga_code import can_we_stop_yet
        okay_to_stop = can_we_stop_yet.evaluate(last_file, second_to_last_file)
        if okay_to_stop:
            pass
        else:
            last_time = os.path.basename(last_file).split("_")[-1]
            new_seed_file = "my_structure_%s" % last_time
            self.clear_ga_vasp_ingredient("vasp_ingredient_match_my_name_in_recipe", new_seed_file)
            self.change_my_status("W")

    def clear_ga_vasp_ingredient(self, ingred_name="", new_seed_file=""):
        """Clear the Genetic Algorithm VASP ingredient.
            Assumes that the VASP ingredient has already
            been evaluated and the result passed
            to the child ingredient.
            Args:
                ingred_name <str>: Ingredient name to clear.
                    Match this carefully with a name from
                    the recipe.
                new_seed_file <str>: New seed file to put
                    in the cleared ingredient.
        """
        if ingred_name == "":
            self.logger.error("Needs an ingredient name!")
            return None
        fullpath = os.path.join(os.path.dirname(self.keywords['name']), ingred_name)
        self.logger.info("Removing directories and files from ingredient specified at %s" % fullpath)
        dircontents = os.listdir(fullpath)
        subfolders = list()
        import shutil
        for diritem in dircontents:
            fulldir = os.path.join(fullpath,diritem)
            if os.path.isdir(fulldir) and diritem.isdigit():
                subfolders.append(fulldir)
        for subfolder in subfolders:
            shutil.rmtree(subfolder)
        files_to_remove = list()
        files_to_remove.append("starting.xyz")
        for filename in files_to_remove:
            os.remove(os.path.join(fullpath, filename))
        shutil.copy(os.path.join(self.keywords['name'], new_seed_file), os.path.join(fullpath, "starting.xyz"))
        from MAST.ingredients import BaseIngredient
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
        structure_list = get_structures(mydir, namelist)
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
            one_structure = pymatgen.io.smartio.read_structure(str_path)
            return_structures.append(one_structure)
        return return_structures
