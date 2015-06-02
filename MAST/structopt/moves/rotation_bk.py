import random
from MAST.structopt.tools.find_defects import find_defects

def rotation(indiv, Optimizer):
    """Move function to perform rotation of a group of atoms
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
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            atms,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
        else:
            atms = indiv[0]
    else:
        atms=indiv[0]
    nat=len(atms)
    if nat != 0:
        if nat<=1:
            natrot2=1
            natrot1=0
        elif nat<=6:
            natrot2=2
            natrot1=0
        else: 
            natrot1=random.randint(1,nat/2)
            natrot2=random.randint(2,(nat-1)/2)
            if natrot2 >= natrot1:
                natrot2 +=1
            else:
                natrot1, natrot2 = natrot2, natrot1
        natrot=natrot2 - natrot1
        atmsr=atms[natrot1:natrot2]
        del atms[natrot1:natrot2]
        ax=['x', '-x','y','-y','z','-z']
        rax=ax[random.randint(0,len(ax)-1)]
        rang=random.uniform(30,180)
        #rang=random.random()*90
        atmsr.rotate(rax,a=rang,center='COM',rotate_cell=False)
        atms.extend(atmsr)
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                atms.extend(indb)
        indiv[0]=atms.copy()
    else:
        natrot=0
        rax=0
        rang=0
    Optimizer.output.write('Rotation Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms rotated = '+repr(natrot)+'\n')
    Optimizer.output.write('Rotation vector = '+repr(rax)+'\n')
    Optimizer.output.write('Rotation angle = '+repr(rang)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='R'+repr(natrot)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv