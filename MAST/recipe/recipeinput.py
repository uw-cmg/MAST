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
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure


ALLOWED_KEYS = {\
                   'recipe_name' : (str, None, 'Recipe Input Name'),\
               }

class RecipeInput(MASTObj):
    """Input to the recipe
    """
    def __init__(self, **kwargs): 
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.recipe_name = self.keywords['recipe_name']
        self.options = InputOptions()
        self.delimiter = ' '
        self.options.set_item("main", "recipe_name", self.recipe_name)

    def addMastInfo(self, program, system_name, scratch_dir=None, section_name="mast"):
        """
           The MAST section input given by the user is added to the
           InputOptions

           Args:
               program : a string used to represent the program
                         to be used by the recipe
               system_name : a string representing the system name

           Returns:
               None

           Raises:

        """
        self.options.set_item(section_name, "program", program)
        self.options.set_item(section_name, "system_name", system_name)

        if scratch_dir is None:
            scratch_dir = os.path.expanduser(os.environ['MAST_SCRATCH'])

        self.options.set_item(section_name, "scratch_directory", scratch_dir)

    def addStructureInfoFromCoords(self, coord_type, coordinates, lattice, addn_info={}, section_name="structure"):
        """
           The Structure section input given by the user is added to the
           InputOptions using the coords mentioned by the user

           Args:
               Coord_type : a string used to represent the coord type
               coordinates : a dict with string of atom_names as key and the 
                              a string for coordinates as value
               lattice : lists of list of strings representing the lattice parameters
               addn_info : a dict of key value pairs used to add to 
                           the input options object
               section_name : usually a constant, which represents the name
                              of the section

           Returns:
               None

           Raises:
        """
        atom_list = list()
        coordinates_list = list()
        for atom, coord in coordinates.iteritems():
            atom_list.append(atom)
            coordinates_list.append(coord.split(self.delimiter))
        np_coordinates = np.array(coordinates_list, dtype="float")

        lattice_list = list()
        for info_str in lattice:
            lattice_list.append(info_str.split(self.delimiter))
        np_lattice = np.array(lattice_list, dtype="float")

        structure = MAST2Structure(lattice=np_lattice,\
                                   coordinates=np_coordinates,\
                                   atom_list=atom_list,\
                                   coord_type=coord_type)
        self.options.set_item(section_name, "structure", structure)
        self.options.set_item(section_name, "coord_type", coord_type)
        self.options.set_item(section_name, "atom_list", atom_list)
        self.options.set_item(section_name, "coordinates", coordinates_list)

        if addn_info:
            for key, value in addn_info.iteritems():
                self.options.set_item(section_name, key, value)
         
    def addStructureInfoFromFile(self, file_name, addn_info={}, section_name="structure"):
        """
           The Structure section input given by the user is added to the
           InputOptions using the coords mentioned by the user

           Args:
               file_name : a string which is a file_name
               addn_info : a dict of key value pairs used to add to 
                           the input options object
               section_name : usually a constant, which represents the name
                              of the section

           Returns:
               None

           Raises:
        """
        if 'poscar' in file_name.lower():
            from pymatgen.io.vaspio import Poscar
            structure = Poscar.from_file(file_name).structure
        elif 'cif' in file_name.lower():
            from pymatgen.io.cifio import CifParser
            structure = CifParser(file_name).get_structures()[0]
        else:
            error = "Cannot build structure from file %s" % file_name
            MASTError(self.__class__.__name__, error)
        
        self.options.set_item(section_name, "structure", structure)

        if addn_info:
            for key, value in addn_info.iteritems():
                self.options.set_item(section_name, key, value)
        
    def addDefectsInfo(self, defects_info, coord_type="cartesian", section_name="defects"):
        """
           Given the input options add the necessary options.

           Args:
               defects_info : a list of defects, each defect item
                              is a list having the type(string), 
                              coordinates(a list of strings) and 
                              symbol(string) is the respective order

               coord_type : a string representing the coordinates type
                            cartesian by default
               section_name : a string representing the section name
                              defects by default
        """
        defect_types = dict()
        defect_types['coord_type'] = coord_type

        count = 0
        for type, coordinates, symbol in defects_info:
            type_dict = defect_types.setdefault('defect%i' % count, dict())

            type_dict['type'] = type
            coordinates = np.array(coordinates, dtype='float')
            type_dict['coordinates'] = coordinates
            type_dict['symbol'] = symbol

            count += 1

        self.options.set_item(section_name, 'num_defects', count)
        self.options.set_item(section_name, 'defects', defect_types)

    def addGlobalIngredientsInfo(self, global_parameters, section_name="ingredients"):
        """
            sets the global parameters for all ingredients

            Args:
                global_parameters : a dict of key value string pairs 
                                    representing the global parameters 
                                    of all ingredients
                section_name : section name which is ingredients by default

            Returns:
                None

            Raises:
        """
        param_dict = global_parameters.copy()
        for key, value in global_parameters.iteritems():
            if key == "mast_kpoints" and type(value) == str:
                kpts = self._get_mast_kpoints(value)
                param_dict[key] = kpts
            else:
                param_dict[key] = value
        self.options.set_item(section_name, "global", global_parameters)

    def _get_mast_kpoints(self, value):
        opt = value.split(self.delimiter)
        kpts = map(int, opt[0].split('x'))
        if len(opt) > 1:
            kpts.append(opt[1])
        return kpts

    def addIngredientsInfo(self, ingredient_name, ingredient_parameters, section_name="ingredients"):
        """
            Adds the ingredients info to the input options

            Args:
                ingredient_name : a string representing the ingredient name
                ingredient_parameters: a dict of key value string pairs which
                                       represents the parameters of the ingredient
                section_name : a string which represents the section name
                               ingredients by default

            Returns:
                None

            Raises:
        """
        global_dict = self.options.get_item(section_name, 'global', dict())
        ingredient_dict = global_dict.copy()

        ingredient_dict.update(ingredient_parameters)

        for key, value in ingredient_dict.iteritems():
            if key == "mast_kpoints" and type(value) == str:
                kpts = self._get_mast_kpoints(value)
                ingredient_dict[key] = kpts

        self.options.set_item(section_name, ingredient_name, ingredient_dict)
            

    def addRecipeInfo(self, recipe_file, addn_info={}, section_name="recipe"): 
        """
           Adds the recipe related options to the input

           Args:
               recipe_file : a string representing the recipe file
               section_name : a string representing the section
                              recipe by default
           
           Returns:
               None

           Raises:
        """
        try:
            recipe_path = os.environ['MAST_RECIPE_PATH']
        except KeyError:
            error = 'MAST_RECIPE_PATH environment variable not set'
            MASTError(self.__class__.__name__, error)

        recipe_file = os.path.join(recipe_path, recipe_file)

        self.options.set_item(section_name, "recipe_file", recipe_file)
        if addn_info:
            for key, value in addn_info.iteritems():
                self.options.set_item(section_name, key, value)


    def addNebInfo(self, hop_info, images=0, section_name="neb"):
        """
           Add Neb related info to the input options

           Args:
               hop_info : a string showing the hop info as below
                          Eg: 1-2 1-3 1-4
               images : an integer representing the number of images
               section_name : a string representing the section name
                              neb by default

            Returns:
                None

            Raises
        """
        hop_list = hop_info.split(self.delimiter)
        hopfrom = dict()

        for hop in hop_list:
            hfrom, hto = [int(hop_elt) for hop_elt in hop.split('-', 1)]
            if not hfrom in hopfrom:
                hopfrom[hfrom] = list()
            hopfrom[hfrom].append(hto)

        self.set_item(section_name, 'images', images)
        self.set_item(section_name, 'hopfrom_dict', hopfrom)


    def addChemicalPotentialsInfo(self, condition_name, condition_info, section_name="chemical_potentials"):
        """
            Adds the chemical potentials info to the input options

            Args:
                condition_name : a string representing the condition name
                                 Eg: As rich
                condition_info : a list of strings representing the conditions
                                 Eg: ['Ga 3.5', 'As 4.5']
                section_name : a string representing the section name
                               chemical potentials by default
        """
        chempot_dict = self.options.get_item(section_name, condition_name, dict())

        cond_dict = dict()
        for condn in condition_info:
            opt = condn.split()
            cond_dict[opt[0]] = float(opt[1])
        chempot_dict[condition_name] = cond_dict

 
