#Fake vasp for testing Amy's GA:
import numpy as np
import time
import pymatgen
import ase
from ase import io as aseio
timestamp = time.asctime()
randenergy = np.random.random_sample()*-100.0 + -200.0
randtime = np.random.random_sample()*100
randpress = np.random.random_sample()*50

#mystruc = pymatgen.io.vaspio.Poscar.from_file("POSCAR").structure
#mystruc.perturb(0.01)
#fake_contcar = pymatgen.io.vaspio.Poscar(mystruc)
#fake_contcar.write_file("CONTCAR",direct=False)
#
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
indiv = aseio.read("POSCAR")
calc = LennardJones()
indiv.set_calculator(calc)
indiv.set_pbc=False
dyn=BFGS(indiv)
try:
    dyn.run(fmax=0.01, steps=1000)
    aseenergy = indiv.get_potential_energy()
    asepress = indiv.get_isotropic_pressure(indiv.get_stress())
except:
    aseenergy = 10000
    asepress = 10000
indiv.set_pbc=True


#poscaratoms.set_cell([4.275, 4.275, 4.275])
ase.io.write("CONTCAR",indiv, "vasp", direct=True, sort=True, vasp5=True)

fake_osz = open("OSZICAR","wb")
fake_osz.write("Output randomly generated at %s\n" % timestamp)
#fake_osz.write("E0=% 3.8f   d E =0.00\n" % randenergy)  
fake_osz.write("E0=% 3.8f   d E =0.00\n" % aseenergy)  
fake_osz.close()

fake_outcar = open("OUTCAR","wb")
fake_outcar.write("Output randomly generated at %s\n" % timestamp)
#fake_outcar.write("external pressure =   %3.8f kB  Pullay stress =  0.00 kB\n" % randpress)
fake_outcar.write("external pressure =   %3.8f kB  Pullay stress =  0.00 kB\n" % asepress)
fake_outcar.write("Pretend we have...reached required accuracy\n")
fake_outcar.write("User time (sec): %3.8f" % randtime)
fake_outcar.close()

