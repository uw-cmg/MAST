from MAST.ingredients.optimize import Optimize
from MAST.utility.metadata import Metadata
from MAST.utility import MASTError
import os
class StaticDefectForNEB(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        label = self.get_my_label()
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname,"parent_structure_" + label)
            self.forward_parent_energy(self.keywords['name'], childname, "parent_energy_" + label)

    def get_my_label(self):
        """Return the defect label, without the word "defect".
            Returns:
                label <str>: defect label
        """
        myname = self.keywords['name']
        mymeta = Metadata(metafile=os.path.join(myname, "metadata.txt"))
        mylabel = mymeta.search_data("defect_label")
        if mylabel == "":
            raise MASTError(self.__class__.__name__,
                "No metadata file for tag defect.")
        plabel = mylabel[1]
        if 'defect_' in plabel:
            plabel = plabel.split("defect_")[-1]
        return plabel
