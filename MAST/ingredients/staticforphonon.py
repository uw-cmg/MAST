from MAST.ingredients.optimize import Optimize

class StaticForPhonon(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        Optimize.update_children(self)
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_extra_restart_files(self.keywords['name'], childname)

