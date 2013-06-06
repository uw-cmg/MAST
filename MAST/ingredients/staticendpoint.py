from MAST.ingredients.optimize import Optimize

class StaticEndpoint(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        mynum = self.get_my_number()
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname,"parent_structure_" + mynum)
            self.forward_parent_energy(self.keywords['name'], childname, "parent_energy_" + mynum)

    def get_my_number(self):
        """For endpoint in the format <sys>_defect<N>_stat, return N.
        """
        tempname = self.keywords['name'].lower()
        numstr = tempname[tempname.find("defect")+6:]
        if numstr.find("defect") == -1:
            pass
        else:
            numstr = numstr[numstr.find("defect")+6:] #allow 'defect' to be in <sys>
        defectstr = numstr[:numstr.find("_")]
        return defectstr

