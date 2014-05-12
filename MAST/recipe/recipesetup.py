##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-05-12 by Zhewen Song
##############################################################
import os, re
import logging
import shlex
from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.utility import loggerutils
from pymatgen.core.structure import Structure
from MAST.recipe import recipeutility as ru
from MAST.recipe.recipeplan import RecipePlan
from MAST.ingredients.baseingredient import BaseIngredient
ALLOWED_KEYS = {
                  'recipeFile'     : (str, None, 'Personalized recipe file'),\
                  'inputOptions'   : (InputOptions, None, 'Input options parsed using input parser'),\
                  'structure'      : (Structure, None, 'Structure to be used to create the ingredient objects'),\
                  'ingredientsDict': (dict, dict(), 'Dictionary of ingredients'),\
                  'workingDirectory': (str, None, 'Working directory')
               }

DATA_PATH = "~/test_dir/"

class RecipeSetup(MASTObj):
    """Parses the personalized recipe file,
        creates the recipe plan object,
        and prepares the ingredient directories.
        Attributes:
            self.recipe_file <str>: Recipe file name
            self.input_options <InputOptions>: Input options for recipe
            self.structure <pymatgen Structure>: Initial structure in recipe
            self.work_dir <str>: Working (recipe) directory)
            self.logger <logging Logger>: MAST monitor and MAST setup-level log
            self.recipe_logger <logging Logger>: Recipe-level log
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.recipe_file    = self.keywords['recipeFile']
        self.input_options  = self.keywords['inputOptions']
        self.structure      = self.keywords['structure']
        self.work_dir    = self.keywords['workingDirectory']
        self.logger = logging.getLogger('mastmon')
        self.logger = loggerutils.add_handler_for_control(self.logger)
        self.recipe_logger = logging.getLogger(self.work_dir)
        self.recipe_logger = loggerutils.add_handler_for_recipe(self.work_dir, self.recipe_logger)

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.work_dir)
        self.recipe_logger.info('Setting up the recipe based on %s' % (self.recipe_file))

    def get_my_ingredient_options(self, name, ingredient_type):
        """Creates the ingredient based on the ingredient type.
            Notes:

            GRJ 6/19/2013: If the ingredient has not be specified in the input
                           file we create it here, and append on the defaults.
                           In addition, the defaults are appended on here, rather
                           than in InputParser (since each ingredient is made here
                           this makes more sense).
            TTM 8/19/13: Unspecified ingredient is "default"
        """
        self.recipe_logger.info('Initializing ingredient %s of type %s' % (name, ingredient_type))
        global_defaults = self.input_options.get_item('ingredients', 'global')

        if (ingredient_type not in self.input_options.get_section_keys('ingredients')):
            self.recipe_logger.info('Ingredient type %s has not be specified in the input file.' % ingredient_type)
            self.recipe_logger.info('Using defaults from ingredients_global.')
            self.input_options.set_item('ingredients', ingredient_type, global_defaults)
        else:
            self.recipe_logger.info('Copying over defaults from ingredients_global for ingredient %s.' % ingredient_type)
            ing_opt = self.input_options.get_item('ingredients', ingredient_type)
            for glob_key, glob_value in global_defaults.items():
                if glob_key not in ing_opt:
                    ing_opt[glob_key] = glob_value
            self.input_options.set_item('ingredients', ingredient_type, ing_opt)

        #self.program = self.input_options.options['ingredients'][ingredient_type]['mast_program']

        ingredient_name = os.path.join(self.work_dir, name)
        pkey_d = self.input_options.get_item('ingredients', ingredient_type).copy()
        mydata = self.metafile.read_data(os.path.basename(ingredient_name)).split(',')
        defect_label=""
        neb_label=""
        charge=""
        phonon_label=""

        if 'defect_' in name.lower():
            for datum in mydata:
                if 'defect_label' in datum:
                    defect_label = datum.split(':')[-1].strip()
                    if not 'defects' in self.input_options.options.keys():
                        raise MASTError(self.__class__.__name__, "No defects section in input file. Error setting up recipe %s." % self.work_dir)
                    defdict = self.input_options.get_item('defects','defects')
                    if not defect_label in defdict.keys():
                        raise MASTError(self.__class__.__name__, "No such label %s found in the defects section dictionary." % defect_label)
                    mydefdict = dict()
                    mydefdict['mast_defect_settings'] = defdict[defect_label]
                    pkey_d.update(mydefdict)
                    break
            if defect_label == "":
                raise MASTError(self.__class__.__name__, "No defect label for %s found in ingredient's metadata.txt" % ingredient_name)
                    
        if 'inducedefect' not in ingredient_type:
            if 'q=' in name.lower():
                for datum in mydata:
                    if 'charge' in datum:
                        charge = int(datum.split(':')[-1])
                        pkey_d['mast_charge'] = charge
                        break

        if 'neb_' in name.lower():
            for datum in mydata:
                if 'neb_label' in datum:
                    neb_label = datum.split(':')[-1].strip()
                    if not 'neb' in self.input_options.options.keys():
                        raise MASTError(self.__class__.__name__, "No neb section in input file. Error setting up recipe %s." % self.work_dir)
                    nebdict = self.input_options.get_item('neb','nebs')
                    if not neb_label in nebdict.keys():
                        raise MASTError(self.__class__.__name__, "No such label %s found in the neb section dictionary." % neb_label)
                    mynebdict = dict()
                    mynebdict['mast_neb_settings'] = nebdict[neb_label]
                    pkey_d.update(mynebdict)
                    break
            if neb_label == "":
                raise MASTError(self.__class__.__name__, "No neb label for %s found in ingredient's metadata.txt" % ingredient_name)

            if 'neb' in name.lower():
                pkey_d.update(self.input_options.options['neb'])
        #if 'phonon' in self.input_options.options.keys():
        #    pkey_d.update(self.input_options.options['phonon'])
        if 'phonon_' in name.lower():
            for datum in mydata:
                if 'phonon_label' in datum:
                    phonon_label = datum.split(':')[-1].strip()
                    def_or_neb_label = phonon_label.split('_')[0]
                    phonon_key = phonon_label.split('_')[-1]
                    defdict = self.input_options.get_item('defects','defects')
                    nebdict = self.input_options.get_item('neb','nebs')
                    phdict=dict()
                    if def_or_neb_label in defdict.keys():
                        phdict['mast_phonon_settings'] = defdict[def_or_neb_label]['phonon'][phonon_key]
                    elif def_or_neb_label in nebdict.keys():
                        phdict['mast_phonon_settings'] = nebdict[def_or_neb_label]['phonon'][phonon_key]
                    else:
                        raise MASTError(self.__class__.__name__, "Neither defect nor neb dictionaries have phonon key %s for ingredient %s." % (phonon_key, name))
                    pkey_d.update(phdict)
                    break
            if phonon_label == "":
                raise MASTError(self.__class__.__name__, "No phonon label for %s found in ingredient's metadata.txt" % ingredient_name)
        if 'chemical_potentials' in self.input_options.options.keys():
            chempotdict=dict()
            chempotdict['mast_chemical_potentials']=self.input_options.options['chemical_potentials']
            pkey_d.update(chempotdict)

        allopt=dict()
        allopt['name'] = ingredient_name
        #allopt['program'] = self.program
        allopt['structure'] = self.structure
        allopt['program_keys'] = pkey_d
        self.recipe_logger.info(pkey_d)
        return allopt
    def get_method_from_ingredient_type(self, ingredtype, methodtype=""):
        """Get write, ready, run, complete, and update methods from
            an ingredient type like volrelax_to_singlerun
            Args:
                ingredtype <str>: ingredient type
                methodtype <str>: method type ("mast_write_method", etc.)
            Returns:
                methodname <str>: method name ("write_default", etc.)
        """
        raise NotImplementedError
        mdict=self.input_options.get_item("ingredients", ingredtype)
        if methodtype in mdict.keys():
            return mdict[methodtype]
        else:
            mdict=self.input_options.get_item("ingredients","global")
            if methodtype in mdict.keys():
                return mdict[methodtype]
            else:
                return None


    def create_recipe_plan(self):
        """Creates a recipe object which has the ingredients and dependency information
        """
        import inspect
        #'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])
        [how_to_update, parents_to_check, how_to_run, recipe_name] = ru.read_recipe(self.recipe_file)

        recipe_obj = RecipePlan(recipe_name, self.work_dir)
        ingredlist = how_to_run.keys()
        for ingred in ingredlist:
            self.update_top_meta_for_ingred(ingred)
            ingredtype = how_to_run[ingred]
            recipe_obj.ingredients[ingred]="I" #set all to initialized
            recipe_obj.ingred_input_options[ingred] = self.get_my_ingredient_options(ingred, ingredtype)
            if not os.path.isdir(os.path.join(self.work_dir, ingred)):
                self.create_ingredient(recipe_obj.ingred_input_options[ingred])
            recipe_obj.write_methods[ingred] = self.get_method_list(ingredtype, "mast_write_method")
            recipe_obj.ready_methods[ingred] = self.get_method_list(ingredtype, "mast_ready_method")
            recipe_obj.run_methods[ingred] = self.get_method_list(ingredtype, "mast_run_method")
            recipe_obj.complete_methods[ingred] = self.get_method_list(ingredtype, "mast_complete_method")
            recipe_obj.update_methods[ingred] = dict()
            for ichild in how_to_update[ingred].keys():
                updingredtype = how_to_update[ingred][ichild]
                recipe_obj.update_methods[ingred][ichild] = self.get_method_list(updingredtype, "mast_update_children_method")
            recipe_obj.parents_to_check = dict(parents_to_check)
        recipe_obj.summary_options = self.get_summary_options()
        return recipe_obj

    def update_top_meta_for_ingred(self, myingred):
        """Update the top metafile for an ingredient.
            Args:
                myingred <str>: ingredient name
        """
        datalist=list()
        datalist.append("ingredient type: %s " % myingred)
        scaling = re.search(r'_\d*x\d*x\d*',myingred)
        if scaling: 
            scalingsize=scaling.group().split('_')[1]
            datalist.append("scaling_size: %s" % scalingsize)
            d_scaling = self.input_options.get_item("structure","scaling")
            kpoints = d_scaling[scalingsize]
            datalist.append("kpoints: %s %s %s" % (kpoints[0],kpoints[1],kpoints[2]) )
        if 'defect_' in myingred:
            if scaling:
                defectlabel = myingred.split('defect_')[1].split('_')[1]
            else: 
                defectlabel = myingred.split('defect_')[1].split('_')[0]
            if defectlabel.isdigit():
                defectlabel = "defect_" + defectlabel
            datalist.append("defect_label: %s" % defectlabel)
        if 'q=' in myingred:
            chargestr = myingred.split('q=')[1].split('_')[0]
            if chargestr[0] == 'p':
                chargelabel=chargestr[1:]
            elif chargestr[0] == 'n':
                chargelabel='-'+chargestr[1:]
            else:
                chargelabel=chargestr
            datalist.append("charge: %s" % chargelabel)
        if 'neb_' in myingred:
            if scaling:
                neblabel = myingred.split('neb_')[1].split('_')[1]
            else:
                neblabel = myingred.split('neb_')[1].split('_')[0]
            datalist.append("neb_label: %s" % neblabel)
        if 'phonon_' in myingred:
            labels = myingred.split('phonon_')[1].split('_')
            try: labels.remove('parse')
            except ValueError: pass
            if scaling: labels.remove(scalingsize)
            phononlabel = '_'.join(labels)
            datalist.append("phonon_label: %s" % phononlabel)
        data=','.join(datalist)
        self.metafile.write_data(myingred, data)

    def create_ingredient(self, my_ingred_input_options):
        """Create the ingredient directory and metadata file.
            Args:
                my_ingred_input_options <Input Options>: single ingredient
                                                            input options
        """
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }


        myingred = BaseIngredient(allowed_keys, 
                        name=my_ingred_input_options['name'],
                        structure=my_ingred_input_options['structure'],
                        program_keys=my_ingred_input_options['program_keys'])
        myingred.write_directory()
        return


    def start(self):
        """Starts the setup process, parse the recipe file
           Use the input options and recipe info to
           create directories and classes required
        """
        if self.recipe_file is None:
            raise MASTError(self.__class__.__name__, "Recipe file not provided!")

        if not os.path.exists(self.recipe_file):
            raise MASTError(self.__class__.__name__, "Recipe file not Found!")

        if not self.input_options:
            raise MASTError(self.__class__.__name__, "Input Options not provided!")


        #print 'DEBUG:, ingredients info =',
        #for ingredient, value in ingredients_info.items():
        #    print ingredient, value
        recipe_plan = self.create_recipe_plan()
        return recipe_plan

    def get_method_list(self, ingredtype, methodtype):
        """Get a method list. Methods should be of
            the form:
                methodname arguments
            where methodname is the method name, and arguments
            are the arguments, as though listed through
            a command line.
            For example, copy_file CONTCAR POSCAR
            Multiple methods may exist with the same keyword
            in the ingredients section, and all will be applied.
            Methods should be separated by a semicolon:
            mast_update_children_method copy_file CONTCAR POSCAR; softlink_file WAVECAR
            Args:
                ingredtype <str>: ingredient type
                methodtype <str>: method type ("mast_write_method", etc.)
            Returns:
                mlist <list>: nested list containing entries
            [[method1name, arg1, arg2,...],[method2name...],...]
        """
        ioptdict=self.input_options.get_item("ingredients", ingredtype)
        if methodtype in ioptdict.keys():
            unparsed = ioptdict[methodtype]
        else:
            globaldict=self.input_options.get_item("ingredients","global")
            if methodtype in globaldict.keys():
                unparsed = globaldict[methodtype]
            else:
                raise MASTError(self.__class__.__name__,"No method type %s in either ingredient type %s or global ingredient." % (methodtype, ingredtype))
        splitmethods = unparsed.split(";")
        mlist = list()
        #self.logger.info(splitmethods)
        for method_arg_item in splitmethods:
            method_arg_item = method_arg_item.strip()
            method_arg_list = shlex.split(method_arg_item)
            methodname = method_arg_list[0].strip()
            arglist = list()
            for argitem in method_arg_list[1:]:
                arglist.append(argitem.strip())
            minputs = list()
            minputs.append(methodname)
            minputs.extend(arglist)
            mlist.append(minputs)
            #if not methodname in mdict.keys():
            #    mdict[methodname]=list(arglist)
            #else:
            #    self.logger.warning("Duplicate entry in method %s for ingredient type %s; ignoring the duplicate." % (methodtype, ingredtype))
        #self.logger.info(mlist)
        return mlist
        #return mdict

    def get_summary_options(self):
        """Get the summary options and give them to the recipe plan."""
        sum_opts = dict()
        if 'summary' in self.input_options.options.keys():
            sumkeys = self.input_options.get_section_keys("summary")
            for sumkey in sumkeys:
                sum_opts[sumkey] = self.input_options.get_item("summary",sumkey)
        return sum_opts
