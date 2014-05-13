import random
from MAST.structopt.tools.find_defects import find_defects

def zp_rotation(indiv, Optimizer):
    """Move function to perform Zero point rotation of atoms
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
    ax=['x', '-x','y','-y','z','-z']
    #Identify random rotation about axis and rotate all atoms
    rax=ax[random.randint(0,len(ax)-1)]
    rang=random.uniform(30,180)
    #rang=random.random()*90
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            atms,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
        else:
            atms = indiv[0]
    else:
        atms=indiv[0]
    atms.rotate(rax,a=rang,center='COM',rotate_cell=False)
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            atms.extend(indb)
    indiv[0]=atms.copy()
    Optimizer.output.write('Zero point full rotation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Rotation vector = '+repr(rax)+'\n')
    Optimizer.output.write('Rotation angle = '+repr(rang)+'\n')
    muttype='ZPR'
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv