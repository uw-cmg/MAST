from MAST.structopt.switches import fitness_switch
import copy

def quench(indiv, Optimizer):
    """Move function to perform a molecular dynamics NVT quenching of structure
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
    try:
        omin = copy.deepcopy(Optimizer.calc.parameters['minimize'])
        mincomd = omin+'\nunfix fix_nve\nfix fix_nvt all nvt temp '+repr(Optimizer.quench_max_temp)+' '+repr(Optimizer.quench_min_temp)+' '+repr(Optimizer.quench_step_size)
        mincomd += '\nrun '+repr(Optimizer.quench_n_steps_1)+'\nunfix fix_nvt\nfix fix_nve all nve\nminimize '+omin
        mincomd += '\nunfix fix_nve\nfix fix_nvt all nvt temp '+repr(Optimizer.quench_max_temp)+' '+repr(Optimizer.quench_min_temp)+' '+repr(Optimizer.quench_step_size)
        mincomd += '\nrun '+repr(Optimizer.quench_n_steps_2)+'\nunfix fix_nvt\nfix fix_nve all nve\nminimize '+omin
    except:
        omin = None
        pms = copy.deepcopy(Optimizer.calc.parameters)
        mincomd = '1.0e-8 1.0e-8 1000 100000\nunfix fix_nve\nfix fix_nvt all nvt temp '+repr(Optimizer.quench_max_temp)+' '+repr(Optimizer.quench_min_temp)+' '+repr(Optimizer.quench_step_size)
        mincomd += '\nrun '+repr(Optimizer.quench_n_steps_1)+'\nunfix fix_nvt\nfix fix_nve all nve\nminimize 1.0e-8 1.0e-8 1000 100000'
        mincomd += '\nunfix fix_nve\nfix fix_nvt all nvt temp '+repr(Optimizer.quench_max_temp)+' '+repr(Optimizer.quench_min_temp)+' '+repr(Optimizer.quench_step_size)
        mincomd += '\nrun '+repr(Optimizer.quench_n_steps_2)+'\nunfix fix_nvt\nfix fix_nve all nve\nminimize 1.0e-8 1.0e-8 1000 100000'
    Optimizer.calc.parameters['minimize'] = mincomd
    Optimizer.calc.parameters['thermosteps'] = Optimizer.lammps_thermo_steps
    outs = fitness_switch([Optimizer,indiv])
    indiv = outs[0]
    Optimizer.output.write('Quench1 NVT Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(outs[1])
    Optimizer.output.write('New fitness : '+repr(indiv.fitness)+'\n')
    muttype='Q'
    if omin == None:
        Optimizer.calc.parameters = pms
    else:
        Optimizer.calc.parameters['minimize'] = omin
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv