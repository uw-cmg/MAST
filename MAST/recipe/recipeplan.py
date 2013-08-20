############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from MAST.ingredients.chopingredient import WriteIngredient
from MAST.ingredients.chopingredient import IsReadyToRunIngredient
from MAST.ingredients.chopingredient import RunIngredient
from MAST.ingredients.chopingredient import IsCompleteIngredient
from MAST.ingredients.chopingredient import UpdateChildrenIngredient
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

    def write_ingredient(self, iname):
        """Write the ingredient files according to the 
            correct method
        """
        methodname = self.write_methods[iname]
        my_ing = WriteIngredient(self.ingred_input_options[iname])
        writeresult=getattr(WriteIngredient, methodname)(my_ing)
        return writeresult

    def complete_ingredient(self, iname):
        """Check if an ingredient is complete
        """
        methodname = self.complete_methods[iname]
        my_ing = IsCompleteIngredient(self.ingred_input_options[iname])
        iscomplete=getattr(IsCompleteIngredient, methodname)(my_ing)
        return iscomplete

    def ready_ingredient(self, iname):
        """Check if an ingredient is ready
        """
        methodname = self.ready_methods[iname]
        my_ing = IsReadyIngredient(self.ingred_input_options[iname])
        isready=getattr(IsReadyIngredient, methodname)(my_ing)
        return isready


    def run_ingredient(self, iname):
        """Run ingredient
        """
        methodname = self.run_methods[iname]
        my_ing = RunIngredient(self.ingred_input_options[iname])
        runresult=getattr(RunIngredient, methodname)(my_ing)
        return runresult

    def update_children(self, iname):
        """Update the children of an ingredient
        """
        upd_results=list()
        for childname in self.update_methods[iname]:
            methodname = self.update_methods[iname][childname]
            my_ing = UpdateChildrenIngredient(self.ingred_input_options[iname])
            updresult=getattr(UpdateChildrenIngredient, methodname)(my_ing)
            upd_results.append(updresult)
        return upd_results

    def check_if_queued_are_complete(self):
        """Check if queued ingredients are complete
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "Q":
                if self.complete_ingredient(iname):
                    self.ingredients[iname] = "C"
                    self.update_children(iname)
        return
    
    def check_if_parents_are_complete(self):
        """Check if parents of waiting ingredients are
            complete.
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "W":
                okay=0
                ptc = list(self.parents_to_check[iname])
                plen = len(ptc)
                for parent in ptc:
                    if self.ingredients[parent] == "C":
                        okay = okay + 1
                if okay == plen:
                    self.ingredients[iname] = "S"

    def run_staged_ingredients(self):
        """Run staged ingredients.
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "S":
                if self.complete_ingredient(iname):
                    self.ingredients[iname] = "C"
                else:
                    if not (self.ready_ingredient(iname)):
                        self.write_ingredient(iname)
                    if self.ready_ingredient(iname):
                        self.run_ingredient(iname)
                        self.ingredients[iname] = "Q"

    def check_recipe_status(self, verbose=1):
        """Check ingredient statuses, and get recipe status
            W = waiting on parents
            I = parents complete, okay to stage
            Q = queued
            C = complete
        Return:
            C = complete
            N = not complete
        """
        total = len(self.ingredients.keys())
        self.check_if_queued_are_complete()
        self.check_if_parents_are_complete()
        self.run_staged_ingredients()
        
        totcomp=0
        totwait=0
        totqueue=0
        totinit=0
        totstage=0
        ilist = self.ingredients.keys()
        ilist.sort()
        if verbose == 1:
            import time
            print time.asctime()
        for iname in ilist:
            if verbose == 1:
                print "%30s = %4s" % (iname, self.ingredients[iname])
            if self.ingredients[iname] == "C":
                totcomp = totcomp + 1
            elif self.ingredients[iname] == "Q":
                totqueue = totqueue + 1
            elif self.ingredients[iname] == "I":
                totinit = totinit + 1
            elif self.ingredients[iname] == "W":
                totwait = totwait + 1
            elif self.ingredients[iname] == "S":
                totstage = totstage + 1
        print "%8s %8s %8s %8s %8s = %8s" % ("INIT","WAITING","STAGED","QUEUED","COMPLETE","TOTAL")
        print "%8i %8i %8i %8i %8i = %8i" % (totinit, totwait, totstage, totqueue, totcomp, total)
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

