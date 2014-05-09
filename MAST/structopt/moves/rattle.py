
def rattle(indiv, Optimizer):
    """Move function to perform rotation of a group of atoms
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
        debug = True/False Boolean
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    indiv[0].rattle(stdev=0.3)
    Optimizer.output.write('Rattle Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='Rat'
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv