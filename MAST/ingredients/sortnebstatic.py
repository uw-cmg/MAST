from MAST.ingredients.neb import NEB
import os 
from MAST.utility.mastfile import MASTFile
from MAST.ingredients.baseingredient import BaseIngredient

class SortNEBStatic(NEB):
    """This class moves static image folders into a new folder. No calculations
        are run.
    """
    def __init__(self, **kwargs):
        NEB.__init__(self, **kwargs)

    def write_files(self):
        NEB.write_files(self)
        dirname = self.keywords['name']
        dirlist = os.listdir(dirname)
        for myfile in dirlist:
           if 'parent_energy' in myfile:
               imfolder = myfile.split('_')[-1]
               if imfolder in dirlist:
                   cpfile = MASTFile(os.path.join(dirname,myfile))
                   cpfile.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, imfolder))


    def is_complete(self):
        dirname = self.keywords['name']
        dirlist = os.listdir(dirname)
        okay=True
        dirct=0
        for mydir in dirlist:
            if os.path.isdir(os.path.join(dirname,mydir)):
                dirct = dirct + 1
                if len(os.listdir(os.path.join(dirname, mydir))) == 2:
                    pass
                else:
                    okay=False
        if dirct == 0:
            return False
        return okay

    def run(self):
        pass
