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
import time

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility.picklemanager import PickleManager

from MAST.recipe.recipeplanner import RecipePlanner

ALLOWED_KEYS = {\
                   'name'              : (str, None, 'Buffet Name'),\
               }


class Buffet(MASTObj):
      """
         Buffet class has a list of dependent recipes and
         relevant feed information between recipes
      """
      def __init__(self, **kwargs):
          MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
          self.name = self.keywords['name']
          self.buffet_plan = list()

      def addRecipe(self, recipe_input, recipe_feeder=None):
          planner = RecipePlanner(recipe_input=recipe_input)
          recipe_plan = planner.plan()
          self.buffet_plan.append((recipe_plan, recipe_feeder.__name__))

      def getBuffetPlan(self):
          plan = []
          for recipe_plan, func_name in self.buffet_plan:
              if func_name is not None:
                  plan.append((recipe_plan, getattr(self, func_name)))
              else:
                  plan.append((recipe_plan, func_name))
          return plan

      def cook(self):
          scratch_dir = os.environ['MAST_SCRATCH']
          pickle_file = os.path.join(scratch_dir, 'buffet_%s_%s.pickle' % (self.name, time.strftime('%Y%m%dT%H%M%S')))
          pm = PickleManager(pickle_file)
          pm.save_variable(self)

      def recipe_feeder(self, completed_recipes, new_recipe):
          """
              User could implement their own feeder or use default ones 
              implemented here.
          """
          pass
