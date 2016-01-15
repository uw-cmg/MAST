import random
import numpy

def atoms_remove(indiv, Optimizer):
    """Move function to alter an individual by removing an atom or group of atoms.
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
   # nato=indiv[0].get_number_of_atoms()-2.0
    nato=50
    atmlist = []
    poplist = []
    atm_pos = indiv[0].get_positions()
    cmax = numpy.maximum.reduce(atm_pos)
    cmin = numpy.minimum.reduce(atm_pos)
    rmax= (cmax-cmin)/2.0
    #zmax = numpy.maximum.reduce(atm_pos)[2] 
    #zmin = numpy.minimum.reduce(atm_pos)[2] 
    for i in range(len(indiv[0])):
       for j in range(3): 
         if atm_pos[i][j] > cmax[j]-0.2*rmax[j] or atm_pos[i][j] < cmin[j] + 0.2*rmax[j] :
            atmlist.append(i)
    list(set(atmlist))
    for i in range(random.randint(1,nato)):
        #indiv[0].pop(random.choice(indiv[0]).index)
        poplist.append(random.choice(atmlist))
    list(set(poplist))
    poplist.sort(reverse=True)
    for i in range(len(poplist)):
        indiv[0].pop(poplist[i])
    Optimizer.output.write('Remove atom mutation\n')
    if Optimizer.structure=='Crystal':
        cell=indiv[0].get_cell()
        a, b, c = cell
        an, bn, cn = [numpy.linalg.norm(v) for v in cell]
        cell[0]=a/an
        cell[1]=b/bn
        cell[2]=c/cn
        indiv[0].set_cell(cell*[(indiv[0].get_number_of_atoms())**0.33*Optimizer.r_ab])
    #muttype='RA'+repr(nato-indiv[0].get_number_of_atoms()+2)
    muttype='RA'+repr(-indiv[0].get_number_of_atoms())
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    Optimizer.output.write(repr(indiv[0].get_positions())+'\n')
    Optimizer.output.write('Removed atoms='+repr(poplist)+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
