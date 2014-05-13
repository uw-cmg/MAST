from MAST.structopt.tools.lammps import LAMMPS
from ase.calculators.vasp import Vasp

def setup_fixed_region_calculator(Optimizer):
    """Function to set up a fixed region calculator for use with LAMMPS and ASE
    """
    calc = self.calc
    pms = copy.deepcopy(calc.parameters)
    nat = sum([c for sym,c,m,mu in self.atomlist])
    try:
        pms['mass'][len(pms['mass'])-1] += '\ngroup RO id >= ' +repr(nat)+'\nfix freeze RO setforce 0.0 0.0 0.0\n'
    except KeyError:
        pms['pair_coeff'][0] += '\ngroup RO id >= '+repr(nat)+'\nfix freeze RO setforce 0.0 0.0 0.0\n'
    ncalc = LAMMPS(parameters=pms, files=calc.files,
    keep_tmp_files=calc.keep_tmp_files, tmp_dir=calc.tmp_dir)
    return ncalc