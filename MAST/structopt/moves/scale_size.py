import random
import math
import numpy

def scale_size(indiv, Optimizer):
    """Move function to change the cell size for crystals or scales size of a cluster
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
    sgn=random.randint(0,1)
    if Optimizer.structure=='Crystal':
        ocell=indiv[0].get_cell()
        cell=ocell
        scale=(sgn+random.random())*ocell
    else:
        ocell=indiv[0].get_cell()
        cell_max=numpy.maximum.reduce(indiv[0].get_positions())
        cell_min=numpy.minimum.reduce(indiv[0].get_positions())
        cell=cell_max-cell_min
        for l in range(len(cell)): 
            cell[l]=math.ceil(cell[l])
            if cell[l]<=1:
                cell[l]=1.5
        indiv[0].set_cell(cell)
        cell=indiv[0].get_cell()
        scale=(sgn+random.random())*cell
    xlen=scale[0][0]-scale[0][1]*scale[1][0]/scale[1][1]
    ylen=scale[1][1]-scale[1][0]*scale[0][1]/scale[0][0]
    zlen=scale[2][2]-scale[2][1]*scale[1][2]/scale[1][1]
    while xlen < 2 or ylen < 2 or zlen < 2:
        sgn=random.randint(0,1)
        scale=(sgn+random.random())*cell
        xlen=scale[0][0]-scale[0][1]*scale[1][0]/scale[1][1]
        ylen=scale[1][1]-scale[1][0]*scale[0][1]/scale[0][0]
        zlen=scale[2][2]-scale[2][1]*scale[1][2]/scale[1][1]
    indiv[0].set_cell(scale, scale_atoms=True)
    if Optimizer.structure != 'Crystal':
        indiv[0].set_cell(ocell)
    Optimizer.output.write('Scaling Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Factor = '+repr(scale)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='SSz'+repr(sgn)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv