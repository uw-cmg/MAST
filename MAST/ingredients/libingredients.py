class BaseIngredient:
    def __init__(self, logger, keywords):
        self.logger    = logger
        self.keywords  = keywords
        self.structure = list() 

    def generate_input(self):
        '''writes the files needed as input for the jobs'''

    def is_complete(self):
        '''Function to check if Ingredient is ready'''


class SinglePoint(BaseIngredient):
    def __init__(self, keywords):
        BaseIngredient.__init__(keywords)



class Optimize(BaseIngredient):
    def __init__(self, keywords):
        BaseIngredient.__init__(keywords)


class InduceDefect(BaseIngredient):
    def __init__(self, keywords):
        BaseIngredient.__init__(keywords)


class OptimizeDefect(BaseIngredient):
    def __init__(self, keywords):
        BaseIngredient.__init__(keywords)

