import pymatgen
import numpy as np
from MAST.ingredients.pmgextend.vasp_extensions import *
mypos = pymatgen.io.vaspio.Poscar.from_file("//home/tam/bin/git/MAST4pymatgen/test/nebtest/ep1/CONTCAR")
mypos.write_file("//home/tam/bin/git/MAST4pymatgen/test/nebtest/POSCAR_no_frozen")
mypos2 = make_one_unfrozen_atom_poscar(mypos, 3)
mypos.write_file("//home/tam/bin/git/MAST4pymatgen/test/nebtest/POSCAR_no_frozen_maybe")
mypos2.write_file("//home/tam/bin/git/MAST4pymatgen/test/nebtest/POSCAR_frozen_3")


