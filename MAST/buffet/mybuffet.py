from MAST.buffet.buffetmanager import Buffet

class MyBuffet(Buffet):
    def __init__(self, name):
        Buffet.__init__(self, name=name)

    def recipe_feeder(self, completed_recipes, new_recipe):
        """
           user defined feeder method which makes the necessary changes to the
           input of the new recipe. The new recipe is a recipe plan which has the
           ingredients objects. The variables of the ingredients could be modified here!!!
        """
        print "recipe feeder not defined !!!"
