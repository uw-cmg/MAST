from MAST.structopt_stem.inp_out.read_xyz import read_xyz
from MAST.structopt_stem.inp_out.write_xyz import write_xyz
from ase import Atom,Atoms
from MAST.structopt_stem.tools.lammps import LAMMPS

# Set up calculator
parcoff = '* * SiC.edip C Si' 
pair_coeff = [parcoff]
mass = ['1 12.011','2 28.0855']
parameters = { 'pair_style' : 'edip', 'pair_coeff' : pair_coeff , 'mass' : mass, 'newton': 'on' }
mincomd = '1e-8 1e-8 5000 10000'
parameters['minimize'] = mincomd
filesL = [ 'SiC.edip' ]
calc = LAMMPS(parameters=parameters, files=filesL)

# Read File
structure = read_xyz('indiv00.xyz',-1)

# Calculate Energy
structure.set_cell([13.092,13.092,13.092])
structure.set_pbc(True)
structure.set_calculator(calc)
OUT=structure.calc.calculate(structure)
totalsol=OUT['atoms']
totalsol.set_pbc(True)
en=OUT['thermo'][-1]['pe']
#Write Relaxed Structure
write_xyz('re-relaxed.xyz',totalsol,repr(en))
fe = en
for sym,u in [('C',-7.371),('Si',-5.3062)]:
	nc=len([atm for atm in structure if atm.symbol==sym])
	fe-= float(nc)*float(u)
print 'SiC Interstitial'
print '    Potential Energy = '+repr(en)
print '    Formation Energy = '+repr(fe)
