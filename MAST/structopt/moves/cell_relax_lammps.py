import os
import copy
from MAST.structopt.tools.lammps import LAMMPS
try:
    from mpi4py import MPI
except ImportError:
    pass

def cell_relax_lammps(indiv, Optimizer):
    """Move function to perform Lammps box/relax for cell. Intended for use in Crystal Optimization
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    cwd1=os.getcwd()
    Optimizer.output.write('LAMMPS fix box/relax performed\n')
    parameters=copy.deepcopy(Optimizer.calc.parameters)
    passflag=True
    try:
        parameters['mass'][len(parameters['mass'])-1] += '\nfix 1 all box/relax iso 0.0 vmax 0.001'
    except KeyError:
        try:
            parameters['pair_coeff'][0] += '\nfix 1 all box/relax iso 0.0 vmax 0.001'
        except KeyError:
            print 'WARNING: LAMMPS Cell relax move trouble with potential. SKIPPING'
            Optimizer.output.write('WARNING: Minimizer quit.  LAMMPS Cell relax move trouble with potential - Box relaxation unsuccessful\n')
            passflag=False
    if passflag:
        filesL = [ Optimizer.pot_file ]
        if Optimizer.lammps_keep_files:
            rank = 0
            path = os.path.join(os.getcwd(),'{0}-rank{1}'.format(Optimizer.filename,rank))
            if not os.path.exists(os.path.join(path,'LAMMPSFiles')):
                            os.mkdir(os.path.join(path,'LAMMPSFiles'))

            real_rank = MPI.COMM_WORLD.Get_rank()
            tmpdir = os.path.join(os.path.join(path, 'LAMMPSFiles'),'rank-{0}'.format(real_rank))
           # calc2 = LAMMPS(parameters=parameters, files=filesL,keep_tmp_files=True, tmp_dir=os.getcwd()+'/'+Optimizer.filename+'/Lammpsrun2/')
            calc2 = LAMMPS(parameters=parameters, files=filesL,keep_tmp_files=True, tmp_dir=tmpdir)
        else:
            calc2 = LAMMPS(parameters=parameters, files=filesL)
        atmsdup=indiv[0].copy()
        atmsdup.set_calculator(calc2)
        try:
            OUT=atmsdup.calc.calculate(atmsdup)
            indiv[0].set_cell(OUT['atoms'].get_cell())
            indiv[0].set_positions(OUT['atoms'].get_positions())
            Optimizer.output.write('Energy = '+repr(OUT['thermo'][-1]['pe']/indiv[0].get_number_of_atoms())+'\n')
        except:
            Optimizer.output.write('WARNING: Minimizer quit - Box relaxation unsuccessful\n')
        #pdb.set_trace()
        os.chdir(cwd1)
        calc2.clean()
        Optimizer.output.write(repr(indiv[0].get_cell())+'\n')
        muttype='LBR'
        if indiv.energy==0:
            indiv.history_index=indiv.history_index+'m'+muttype
        else:
            indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
