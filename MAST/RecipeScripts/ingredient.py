class Ingredient:
    def __init__(self, name, ingredient_type, input_options):
        self.name            = name
        self.ingredient_type = ingredient_type
        self.obj             = self.create_ingredient(input_options)

    def create_ingredient(self, input_options):
        '''Creates the ingredient based on the ingredient type'''

    def setup(self):
        '''Create the directory and call the ingredient's
           generate input function
        '''


