############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from MAST.ingredients.chopingredient import ChopIngredient
class RecipePlan:
    """Contains the entire recipe plan. It consists of
       - Ingredients Objects
       - Dependency Dict.
    """
    def __init__(self, name):
        self.name            = name
        self.ingredients     = dict()  #name=status
        self.update_methods   = dict()
        self.parents_to_check= dict()
        self.run_methods      = dict()
        self.write_methods    = dict()
        self.ready_methods    = dict()
        self.complete_methods = dict()
        self.ingred_input_options = dict()
        #print 'GRJ DEBUG: Initializing RecipePlan'

    def write_ingredient(iname):
        """Write the ingredient files according to the 
            correct method
        """
        methodname = self.write_methods[iname]
        my_ing = ChopIngredient(self.ingred_input_options[iname])
        writeresult=getattr(ChopIngredient, methodname)(my_ing)
        return writeresult

    def complete_ingredient(iname):
        """Check if an ingredient is complete
        """
        methodname = self.complete_methods[iname]
        my_ing = ChopIngredient(self.ingred_input_options[iname])
        iscomplete=getattr(ChopIngredient, methodname)(my_ing)
        return iscomplete

    def ready_ingredient(iname):
        """Check if an ingredient is ready
        """
        methodname = self.ready_methods[iname]
        my_ing = ChopIngredient(self.ingred_input_options[iname])
        isready=getattr(ChopIngredient, methodname)(my_ing)
        return isready


    def run_ingredient(iname):
        """Run ingredient
        """
        methodname = self.run_methods[iname]
        my_ing = ChopIngredient(self.ingred_input_options[iname])
        runresult=getattr(ChopIngredient, methodname)(my_ing)
        return runresult

    def update_children(iname):
        """Update the children of an ingredient
        """
        upd_results=list()
        for childname in self.update_methods[iname]:
            methodname = self.update_methods[iname][childname]
            my_ing = ChopIngredient(self.ingred_input_options[iname])
            updresult=getattr(ChopIngredient, methodname)(my_ing)
            upd_results.append(updresult)
        return upd_results

    def check_if_queued_are_complete():
        """Check if queued ingredients are complete
        """
        for iname in self.ingredients.keys():
            if self.ingredient[iname] == "Q":
                if self.complete_ingredient(iname):
                    self.ingredient[iname] = "C"
                    self.update_children(iname)
        return
    
    def check_if_parents_are_complete():
        """Check if parents of waiting ingredients are
            complete.
        """
        for iname in self.ingredients.keys():
            if self.ingredient[iname] == "W":
                okay=0
                ptc = list(self.parents_to_check[iname])
                plen = len(ptc)
                for parent in ptc:
                    if self.ingredients[parent] == "C":
                        okay = okay + 1
                if okay == plen:
                    self.ingredients[iname] = "S"

    def run_staged_ingredients():
        """Run staged ingredients.
        """
        for iname in self.ingredients.keys():
            if self.ingredient[iname] == "S":
                self.write_ingredient(iname)
                if self.ready_ingredient(iname):
                    self.run_ingredient(iname)
                    self.ingredient[iname] = "Q"


    def check_recipe_status():
        """Check ingredient statuses, and get recipe status
            W = waiting on parents
            I = parents complete, okay to stage
            Q = queued
            C = complete
        Return:
            C = complete
            R = running
        """
        total = len(self.ingredients.keys())
        self.check_if_queued_are_complete()
        self.check_if_parents_are_complete()
        self.run_staged_ingredients()
        
        totcomp=0
        for iname in self.ingredients():
            if self.ingredients[iname] == "C":
                totcomp = totcomp + 1
        if totcomp == total:
            return True
        return False

    def add_ingredient(self, ingredient_name, ingredient):
        """Used to add an ingredient_object corresponding to an ingredient name
        """
        self.ingredients[ingredient_name] = ingredient

    def get_ingredient(self, ingredient_name):
        """Used to get an ingredient_object corresponding to an ingredient name
        """
        return self.ingredients.get(ingredient_name)

    #def add_parent(self, ingredient_name, parent_name):
    #    """Used to build the dependency dict using the parent child relationship
    #    """
    #    self.dependency_dict.setdefault(ingredient_name, list()).append(parent_name)

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

