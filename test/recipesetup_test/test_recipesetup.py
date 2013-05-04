import unittest
import os

from MAST.utility import InputOptions
from MAST.utility import MASTError
from MAST.recipe.recipesetup import RecipeSetup
from MAST.ingredients.ingredients_loader import IngredientsLoader

class TestRecipeSetup(unittest.TestCase):
    def setUp(self):
        if not os.path.exists("test_directories"):
            os.mkdir("test_directories")
        self.input_options = InputOptions()
        self.input_options.set_item('mast', 'program', 'vasp')
        self.input_options.set_item('mast', 'working_directory', 'test_directories/')

        self.ing_loader = IngredientsLoader()
        self.ing_loader.load_ingredients()
        self.ingredients_dict = self.ing_loader.ingredients_dict

    def tearDown(self):
        try:
            for dir in os.listdir('test_directories/'):
                os.rmdir('test_directories/%s' % dir)
        except OSError:
            pass

    def test_unprovided_recipe(self):
        setup_obj = RecipeSetup(recipeFile=None, inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        self.assertRaises(MASTError, setup_obj.start)

    def test_unfound_recipe(self):
        setup_obj = RecipeSetup(recipeFile="dummy_file", inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        self.assertRaises(MASTError, setup_obj.start)

    def test_unprovided_input_options(self):
        setup_obj = RecipeSetup(recipeFile="test_personal_recipe.txt", inputOptions=None, structure=None, ingredientsDict=self.ingredients_dict)
        self.assertRaises(MASTError, setup_obj.start)

    def test_unprovided_ingredients_dict(self):
        setup_obj = RecipeSetup(recipeFile="test_personal_recipe.txt", inputOptions=self.input_options, structure=None, ingredientsDict=None)
        self.assertRaises(MASTError, setup_obj.start)

    def test_unimplemented_ingredient(self):
        setup_obj = RecipeSetup(recipeFile="unimplemented_recipe.txt", inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        self.assertRaises(MASTError, setup_obj.start)

    def test_undefined_ingredient(self):
        setup_obj = RecipeSetup(recipeFile="undefined_recipe.txt", inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        self.assertRaises(MASTError, setup_obj.start)

    def test_parse_recipe(self):
        setup_obj = RecipeSetup(recipeFile="test_personal_recipe.txt", inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        ingredients_info, recipe_name = setup_obj.parse_recipe()
        self.assertEqual(ingredients_info, {'testSys_perfect_opt1': ['optimize', {'testSys_perfect_opt2': ['structure', 'neb']}], 'testSys_perfect_opt2': ['optimize', {}]})

    def test_start(self):
        setup_obj = RecipeSetup(recipeFile="test_personal_recipe.txt", inputOptions=self.input_options, structure=None, ingredientsDict=self.ingredients_dict)
        recipe_plan = setup_obj.start()
        self.assertTrue(recipe_plan != None)
        self.assertEqual(recipe_plan.name, "test") 
        self.assertEqual(recipe_plan.dependency_dict, {'testSys_perfect_opt2' : ['testSys_perfect_opt1']})
        expected_dirs = ['testSys_perfect_opt1', 'testSys_perfect_opt2']
        for dir in os.listdir("test_directories/"):
            self.assertTrue(dir in expected_dirs)
