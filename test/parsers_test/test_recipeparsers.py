import unittest
import os
import filecmp

from MAST.utility import InputOptions
from MAST.utility import MASTError
from MAST.parsers.recipeparsers import RecipeParser 

class TestRecipeParser(unittest.TestCase):
    def setUp(self):
        self.input_options = InputOptions()
        self.input_options.set_item("mast", "system_name", "testSys")
        self.input_options.set_item("defects", "num_defects", 2)
        self.input_options.set_item("neb", "images", 2)
        self.input_options.get_item("neb", "hopfrom_dict", {1:[2]})

    def tearDown(self):
        try:
            os.remove('test_simple_personal_recipe.txt')
            os.remove('test_complex_personal_recipe.txt')
        except OSError:
            pass
 
    def test_unprovided_template(self):
        recipe_obj = RecipeParser(inputOptions=self.input_options, personalRecipe='test_simple_personal_recipe.txt')
        self.assertRaises(MASTError, recipe_obj.parse)  

    def test_unfound_template(self):
        recipe_obj = RecipeParser(templateFile='unfound_test_recipe.txt', inputOptions=self.input_options)
        self.assertRaises(MASTError, recipe_obj.parse)

    def test_unprovided_personal_recipe(self):
        recipe_obj = RecipeParser(templateFile='test_recipe_simple.txt', inputOptions=self.input_options)
        self.assertRaises(MASTError, recipe_obj.parse)

    def test_unprovided_input_options(self):
        recipe_obj = RecipeParser(templateFile='test_recipe_simple.txt', personalRecipe='test_simple_personal_recipe.txt')
        self.assertRaises(MASTError, recipe_obj.parse)

    def test_parse_simple_run(self):
        recipe_obj = RecipeParser(templateFile='test_recipe_simple.txt', inputOptions=self.input_options, personalRecipe='test_simple_personal_recipe.txt')
        recipe_obj.parse()
        self.assertTrue(filecmp.cmp('test_simple_personal_recipe.txt', 'simple_personal_recipe.txt')) 

    def test_parse_complex_run(self):
        recipe_obj = RecipeParser(templateFile='test_recipe_complex.txt', inputOptions=self.input_options, personalRecipe='test_complex_personal_recipe.txt')
        recipe_obj.parse()
        self.assertTrue(filecmp.cmp('test_complex_personal_recipe.txt', 'complex_personal_recipe.txt')) 
