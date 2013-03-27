import pymatgen
import os
from MAST.ingredients.performneb import PerformNEB
os.chdir("//home/tam/bin/git/MAST4pymatgen/test/nebtest")
NEBing=PerformNEB(name="nebtest_neb10-11", program="no program", program_keys={'images':3})
NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
image_list=NEBing.write_files()

print "Complete? :", NEBing.is_complete()
import shutil
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/01")
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/02")
print "Complete? :", NEBing.is_complete()
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/03")
print "Complete? :", NEBing.is_complete()

