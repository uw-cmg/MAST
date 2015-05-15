import os
from ase import Atom, Atoms
from ase.io import read
from MAST.structopt.inp_out import read_xyz, write_xyz

def gen_solid(solidfile,solidcell,outfilename,calc=False,calcmeth=None):
    """Function to load a bulk solid from a file for use in Defect structure optimization
    Inputs:
        solidfile=String of filename to load
        solidcell=List/Matrix of cell parameters for ASE Atoms class
        outfilename=String of filename to write solid
        calc=False/calculator object for evaluating energy of solid
        calcmeth='VASP' or other method for calculating the energy of the solid
    Outputs:
        solid as ASE Atoms class
        energy and string if calc is not false
    """
    try:
        sol = read_xyz(solidfile)
    except Exception as e1:
        try:
            sol = read(solidfile)
        except Exception as e2:
            raise RuntimeError('Encountered errror:'+repr(e1)+' '+repr(e2)+
                ' While trying to read solid file given as:'+repr(solidfile))
    cell = solidcell
    sol.set_cell(cell)
    sol.set_pbc(True)
    #Evaluate pure Bulk structure
    if calc:
        cwd = os.getcwd()
        sol.set_calculator(calc)
        stro = ''
        try:
            if calcmeth == 'VASP':
                en = sol.get_potential_energy()
                calcb = Vasp(restart=True)
                sol = calcb.get_atoms()
                PureBulkEnpa = en/sol.get_number_of_atoms()
            else:
                OUT = sol.calc.calculate(sol)
                PureBulkEnpa = OUT['thermo'][-1]['pe']/sol.get_number_of_atoms()
                sol = OUT['atoms']
                sol.set_pbc(True)
        except:
            stro = 'WARNING: Unable to calculate energy of pure bulk solid'
            PureBulkEnpa = 0
        os.chdir(cwd)
        # Write bulk file to directory
        write_xyz(outfilename,sol,PureBulkEnpa)
        return sol, PureBulkEnpa, stro
    else:
        # Write bulk file to directory
        write_xyz(outfilename,sol,'Pure Bulk')
        return sol