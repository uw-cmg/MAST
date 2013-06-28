from MAST.ingredients.optimize import Optimize
from MAST.utility import MASTError

class StaticDefect(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)

    def update_children(self):
        mynum = self.get_my_label()
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname,"parent_structure_" + mynum)
            self.forward_parent_energy(self.keywords['name'], childname, "parent_energy_" + mynum)

    def get_my_label(self, signal="defect", altsignal="image"):
        """For defect in the format 
            <sys>_..._<signal>_<label>_..._stat, return label.
            Args:
                signal <str>: Signal tag that precedes the label and is
                                separated from the label by an underscore.
                                Only the FIRST OCCURRENCE of the signal is
                                noticed.
                altsignal <str>: Alternative signal tag 
        """
        import os
        tempname = os.path.basename(self.keywords['name'].lower())
        tempsplit = tempname.split('_')
        if signal in tempsplit:
            defectidx = tempsplit.index(signal)
            try:
                defectstr = tempsplit[defectidx+1]
            except IndexError:
                raise MASTError(self.__class__.__name__,
                    "Defect signal '%s_' is not followed by a label" % signal)
        elif altsignal in tempsplit:
            defectidx = tempsplit.index(altsignal)
            defectstr = '_'.join(tempsplit[defectidx+1:defectidx+3])
            if defectstr == "":
                raise MASTError(self.__class__.__name__,
                    "Defect signal '%s_' is not followed by a label" % altsignal)
        else:
            raise MASTError(self.__class__.__name__,
                "Defect label in ingredient name is not preceded by '%s_' or '%s_' " % (signal, altsignal))
        return defectstr

