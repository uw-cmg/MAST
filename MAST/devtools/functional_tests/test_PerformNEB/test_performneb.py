import pymatgen
import os
#os.chdir("//home/tam/bin/git/tam_k3/MAST/ingredients")
from MAST.MAST.ingredients.performneb import PerformNEB
os.chdir("//home/tam/bin/git/MAST4pymatgen/devtools/functional_tests/test_PerformNEB")
#contcar_initial=pymatgen.io.vaspio.Poscar.from_file("//home/tam/test_PerformNEB/ep1FCC",False)
#contcar_final=pymatgen.io.vaspio.Poscar.from_file("//home/tam/test_PerformNEB/ep2FCC",False)
#NEBing=PerformNEB(dir_name="my_test_neb", structure_init=contcar_initial.structure, structure_final=contcar_final.structure, images=4)
NEBing=PerformNEB(dir_name="my_test_neb", parent_init="ep1", parent_final="ep2", program="no program", images=4)
NEBing=PerformNEB(dir_name="my_test_neb", parent_init="ep1", parent_final="ep2", program="vasp", images=4)
image_list=NEBing.generate_files()

print "Complete? :", NEBing.is_complete()
import shutil
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/01")
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/02")
print "Complete? :", NEBing.is_complete()
shutil.copy(os.getcwd() + "/OUTCAR", "my_test_neb/03")
print "Complete? :", NEBing.is_complete()

