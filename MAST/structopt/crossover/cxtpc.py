import random
import copy

def cxtpc(ind1, ind2, Optimizer):
    """ Simple Two point Crossover for crystals
    Positions exchanged as fraction of cell
    Maintains concentration
    Maintains total number of atoms
    Random clusters not necessarily geometry based
    """
    if 'CX' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Two Pt Crystal-Fraction Cx between individual '
        +repr(ind1.index)+' and individual '+repr(ind2.index)+'\n')
    indi1 = ind1[0]
    indi2 = ind2[0]
    pos1 = indi1.get_scaled_positions()
    pos2 = indi2.get_scaled_positions()
    size = min(len(pos1), len(pos2))
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else:
        # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1
    posi1 = copy.deepcopy(pos1)
    posi2 = copy.deepcopy(pos2)
    posi1[cxpoint1:cxpoint2] = pos2[cxpoint1:cxpoint2]
    posi2[cxpoint1:cxpoint2] = pos1[cxpoint1:cxpoint2]
    indi1.set_scaled_positions(posi1)
    indi2.set_scaled_positions(posi2)
    if debug:
        print 'DEBUG CX: CXPT1 = ',cxpoint1,'CXPT2 = ',cxpoint2
    ind1[0] = indi1
    ind2[0] = indi2
    return ind1, ind2

