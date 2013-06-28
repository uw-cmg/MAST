############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################

class RecipePlan:
    """Contains the entire recipe plan. It consists of
       - Ingredients Objects
       - Dependency Dict.
    """
    def __init__(self, name):
        self.name            = name
        self.ingredients     = dict()
        self.dependency_dict = dict()
        #print 'GRJ DEBUG: Initializing RecipePlan'

    def add_ingredient(self, ingredient_name, ingredient):
        """Used to add an ingredient_object corresponding to an ingredient name
        """
        self.ingredients[ingredient_name] = ingredient

    def get_ingredient(self, ingredient_name):
        """Used to get an ingredient_object corresponding to an ingredient name
        """
        return self.ingredients.get(ingredient_name)

    def add_parent(self, ingredient_name, parent_name):
        """Used to build the dependency dict using the parent child relationship
        """
        self.dependency_dict.setdefault(ingredient_name, list()).append(parent_name)

    def ingredient_iterator(self):
        """Iterates through the ingredients dict and returns the ingredients one by one
        """
        for ingredient_name, ingredient_obj in self.ingredients.iteritems():
            yield ingredient_name, ingredient_obj

    def __iter__(self):
        """Iterates through the ingredients dict and returns the ingredients one by one
        """
        for ingredient_name, ingredient_obj in self.ingredients.iteritems():
            yield ingredient_obj

