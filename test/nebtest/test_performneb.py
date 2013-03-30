import pymatgen
import os
from MAST.ingredients.performneb import PerformNEB
os.chdir("//home/tam/bin/git/MAST4pymatgen/test/nebtest")
NEBing=PerformNEB(dir_name="my_test_neb", parent_init="ep1", parent_final="ep2", program="no program", images=3)
NEBing=PerformNEB(dir_name="my_test_neb", parent_init="ep1", parent_final="ep2", program="vasp", images=3)
image_list=NEBing.generate_files()

print "Complete? :", NEBing.is_complete()
import shutil
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/01")
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/02")
print "Complete? :", NEBing.is_complete()
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/03")
print "Complete? :", NEBing.is_complete()

