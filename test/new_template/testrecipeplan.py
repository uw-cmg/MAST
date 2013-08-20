from MAST.recipe import recipeutility as ru
from MAST.recipe.recipesetup import RecipeSetup
from MAST.recipe.recipeplan import RecipePlan
from MAST.parsers.inputparser import InputParser
from MAST.mast import MAST
#ru.read_recipe("//home/tam/tammast/recipe_templates/treetest.txt")

mast_obj = MAST(inputfile="//home/tam/tammast/test/new_template/test_small_new_template.inp")
mast_obj.check_independent_loops()
rs = RecipeSetup(recipeFile = "//home/tam/tammast/recipe_templates/treetest.txt", inputOptions = mast_obj.input_options,structure=mast_obj.input_options.get_item('structure','structure'))
rplan = rs.start()
rplan.check_recipe_status()

