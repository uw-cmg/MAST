from MAST.ingredients.neb import NEB

class NEBLowMesh(NEB):
    def __init__(self, **kwargs):
        PerformNEB.__init__(self, **kwargs)
    
    def update_children(self):
        """Update to ANOTHER NEB."""
        import os
        for childname in self.keywords['child_dict'].iterkeys():
            myct=1
            while myct <= self.keywords['program_keys']['images']:
                imno = str(myct).zfill(2)
                impath = os.path.join(self.keywords['name'], imno)
                self.forward_parent_structure(impath, childname,"parent_structure_" + '-'.join(self.labels) + '_' + imno)
                myct = myct + 1
