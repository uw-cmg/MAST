import os
import random
from ase import Atom, Atoms
from ase.optimize import BFGS
from MAST.structopt.tools.lammps import LAMMPS
from MAST.structopt.tools.setup_calculator import setup_calculator

def ase_minimization(indiv, Optimizer):
    """Function to use built in ASE minimizers to minimize atomic positions in structure.
    Input:
        indiv = Individual class object to be optimized
        Optimizer = Optimizer class object with needed parameters
    Output:
        indiv = Optimized Individual class object.
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    cwd1=os.getcwd()
    olammpsmin = Optimizer.lammps_min
    if Optimizer.lammps_min:
        Optimizer.lammps_min = None
    calc2 = setup_calculator(Optimizer)
#     if 'mass' in indiv[0].get_calculator():
#         mass2 = ['1 '+ str(Optimizer.atomlist[0][2])]
#         if len(Optimizer.atomlist) > 1:
#             for i in range(len(Optimizer.atomlist)-1):
#                 mass2.append(str(i+2) + ' ' + str(Optimizer.atomlist[i+1][2]))
#         calc2=LAMMPS(parameters={ 'pair_style' : Optimizer.pair_style, 'pair_coeff' : Optimizer.pair_coeff , 'mass' : mass2 },files=[ Optimizer.pot_file ])
#     else:
#         calc2=LAMMPS(parameters={ 'pair_style' : Optimizer.pair_style, 'pair_coeff' : Optimizer.pair_coeff},files=[ Optimizer.pot_file ])
    if Optimizer.structure==Defect:
        nat=indiv[0].get_number_of_atoms
        sol=Atoms()
        sol.extend(indiv[0])
        sol.extend(indiv.bulko)
        sol.set_calculator(calc2)
        sol.set_cell(indiv.bulko.get_cell())
        sol.set_pbc(True)
        dyn=BFGS(sol)
        dyn.run(fmax=0.001, steps=2500)
        positions=sol[0:nat].get_positions()
        indiv[0].set_positions(positions)
    else:
        atomsdup=indiv[0].copy()
        atomsdup.set_calculator(calc2)
        dyn=BFGS(indiv[0])
        dyn.run(fmax=0.001, steps=2500)
        positions=atomsdup.get_positions()
        indiv[0].set_positions(positions)
    os.chdir(cwd1)
    calc2.clean()
    Optimizer.lammps_min = olammpsmin
    Optimizer.output.write('ASE Minimization mutation performed on individual = '+repr(indiv.index)+'\n')
    muttype='ASEM'
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
		
		