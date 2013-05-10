from MAST.utility import MASTError

class RecipeBatch:
    """The recipe batch contains a list of dependent recipes, in order.
        (Independent recipes should each have their own recipe batch.)
    """
    def __init__(self, name):
        self.recipe_plan_list = list() #list of recipe plans
        self.recipe_plan_results = list()
        self.name = name #This could also be an ID number??
        self.current_recipe_plan_index = 0 #index of current recipe plan
        self.batch_complete = False

    def update_recipe_plan(self):
        """Set self.current_recipe_plan_index to the current recipe,
            which is the first incomplete recipe in the list.
            ***Is the recipe plan a pickle or not a pickle?
            Also set self.batch_complete if the recipe batch is complete.
        """
        for recipe_plan in recipe_list:
            for ingredient in recipe_plan.ingredient_iterator():
                if not ingredient.is_complete():
                    self.current_recipe_plan_index = self.recipe_plan_list.index(recipe_plan)
                    return
        self.current_recipe_plan_index = len(self.recipe_plan_list)
        self.batch_complete = True


    def save_recipe_plan_output(self, recipe_plan):
        """Save the output for one recipe plan into self.previous_results.
            This output could be a dictionary with entries for structures,
            energies, etc.
        """
        return

    def get_relevant_complete_ingredient(self, recipe_plan):
        """Get the relevant complete ingredient for a completed recipe plan.
            This is often but not always the last ingredient.
        """
        return complete_ingredient

    def get_ingredient_results(self, complete_ingredient):
        """Save the results from a relevant complete ingredient which are
            necessary for the next recipe.
        """
        myresults = dict()
        myresults['structure'] = Structure()
        myresults['energy'] = complete_ingredient.get_my_energy() # use pymatgen here...
        return myresults

    def save_recipe_plan_results(self, recipe_plan):
        """Save the results from a completed recipe plan's relevant complete
            ingredient.
        """
        complete_ingredient = self.get_relevant_complete_ingredient(recipe_plan)
        complete_results = self.get_ingredient_results(complete_ingredients)
        self.recipe_plan_results.append(complete_results)
        return

    def get_ingredient_for_update(self, recipe_plan):
        """Get the ingredient to be updated in the next recipe plan.
            This should really be the first ingredient...
        """
        return first_ingredient

    def compare_two_recipe_plans(self, recipe_plan_1, recipe_plan_2):
        """Compare two dependent recipes to see if some ending criterion 
            has been met, for example, the energy of the last complete
            recipe and the second to last complete recipe. 
        """
        return

    def update_ingredient(self, next_ingredient):
        """Update the correct ingredient in the next recipe plan with
            the correct information.
        """
        #to update a structure, something like:
        newstructure = self.recipe_plan_results[self.current_recipe_plan_index - 1]['structure']
        forward_structure(next_ingredient, structure)
        #we may instead of this kind of updating want recipes to create
        #subsequent INPUT FILES and then call mast, in order to create a new
        #recipe
        return

    def forward_structure(self, next_ingredient, structure):
        """Forward a structure to the next ingredient.
        """
        program = next_ingredient.keywords['program']
        if program == 'vasp':
            Poscar(structure)
            Poscar.write_file(os.path.join(next_ingredient.keywords['name'],'POSCAR'))
            #create a poscar with this structure
        #we may instead of this kind of updating want recipes to create
        #subsequent INPUT FILES and then call mast, in order to create a new
        #recipe
        return
