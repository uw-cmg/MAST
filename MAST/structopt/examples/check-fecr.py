from MAST.structopt.io.read_xyz import read_xyz
from MAST.structopt.io.write_xyz import write_xyz
from ase import Atom,Atoms
from MAST.structopt.tools.lammps import LAMMPS
import os

# Set up calculator
parcoff = '* * FeCr.cdeam Cr Fe' 
pair_coeff = [parcoff]
parameters = { 'pair_style' : 'eam/cd', 'pair_coeff' : pair_coeff }
mincomd = '1e-8 1e-8 5000 10000'
parameters['minimize'] = mincomd
filesL = [ 'FeCr.cdeam' ]
calc = LAMMPS(parameters=parameters, files=filesL, keep_tmp_files=True, tmp_dir=os.path.join(os.getcwd(), 'LAMMPSFiles'))

files = os.listdir(os.getcwd())
files = [file for file in files if '.xyz' in file]
for file in flist:
    # Read File
    structure = read_xyz(file,-1)
    # Calculate Energy
    structure.set_cell([22.96, 22.96, 22.96])
    structure.set_pbc(True)
    structure.set_calculator(calc)
    OUT=structure.calc.calculate(structure)
    totalsol=OUT['atoms']
    totalsol.set_pbc(True)
    en=OUT['thermo'][-1]['pe']
    #Write Relaxed Structure
    write_xyz('r{0}'.format(file),totalsol,repr(en))
    fe = en
    for sym,u in [('Cr',-3.8363),('Fe',-4.21224)]:
        nc=len([atm for atm in structure if atm.symbol==sym])
        fe-= float(nc)*float(u)
    print 'FeCr Interstitial - ' + file
    print '    Potential Energy = '+repr(en)
    print '    Formation Energy = '+repr(fe)