##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import shutil
import sys
import time
from MAST.ingredients.baseingredient import BaseIngredient
from pymatgen.core.structure import Structure
class CustomChopIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        """Please not modify this init method."""
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)


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


