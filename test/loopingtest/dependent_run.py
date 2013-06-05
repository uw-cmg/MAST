from MAST.recipe.recipeinput import RecipeInput
from MAST.buffet.mybuffet import MyBuffet

####################################
# input parameters
####################################

#mast info
PROGRAM     = "VASP"
SYSTEM_NAME = "MgAl"

#structure info
COORD_TYPE  = "fractional"
COORDINATES = {\
                'Mg' : '0.000000 0.000000 0.000000',\
                'Al' : '0.500000 0.500000 0.500000',\
              }
              
LATTICE     = [\
                 "3.0 0.0 0.0",\
                 "0.0 3.0 0.0",\
                 "0.0 0.0 3.0"\
              ]

#ingredients info
MAST_KPOINTS = ['2x2x2', '3x3x3']
ING_GLOBAL_PARAM = {\
                       'mast_kpoints' : None,\
                       'xc'           : 'pbe',\
                   }

#optimize ingredient info
OPTIMIZE_ING  = "optimize"
OPTIMIZE_INFO = {\
                    'encut'  : '300',\
                    'ibrion' : '2',\
                }

#recipe info
RECIPE_FILE = 'recipe_test.txt'

####################################


buffet_obj = MyBuffet(name="dependent_recipe")

#dependent run
for run_idx in xrange(len(MAST_KPOINTS)):
    recipe_input = RecipeInput(recipe_name="dependent_recipe_%s" % MAST_KPOINTS[run_idx])

    #add input parameters
    recipe_input.addMastInfo(PROGRAM, SYSTEM_NAME)
    recipe_input.addStructureInfoFromCoords(COORD_TYPE, COORDINATES, LATTICE)

    #change the required input parameter
    ING_GLOBAL_PARAM['mast_kpoints'] = MAST_KPOINTS[run_idx]

    recipe_input.addGlobalIngredientsInfo(ING_GLOBAL_PARAM)
    recipe_input.addIngredientsInfo(OPTIMIZE_ING, OPTIMIZE_INFO)
    recipe_input.addRecipeInfo(RECIPE_FILE)

    #add recipe input to buffet and cook it
    buffet_obj.addRecipe(recipe_input, buffet_obj.recipe_feeder)

buffet_obj.cook()


