from ase import Atoms, Atom
from MAST.structopt.inp_out.write_xyz import write_xyz

def write_pop(Optimizer,pop):
    """
    Function to write the data from a given population.
    Inputs:
        Optimizer = structopt Optimizer class object with output file data
        pop = List of structopt Individual class objects with data to be written
    Outputs:
        None. Data is written to output files provided by Optimizer class object
    """
    Optimizer.output.write('\n--New Population--\n')
    #Write structures
    for ind in pop:
        update_outfile(ind, Optimizer.output)
        if Optimizer.genealogy:
            update_genealogy(ind, Optimizer.Genealogyfile)
        if Optimizer.swaplist: 
            swaplist_check(ind,Optimizer.structure,Optimizer.output)
        if Optimizer.indiv_defect_write:
            write_xyz(Optimizer.ifiles[ind.index],ind[0],ind.energy)
        update_structsumfile(ind, Optimizer.files[Optimizer.nindiv])
        positions = update_structfile(ind, Optimizer.files[ind.index], Optimizer)
        Optimizer.output.write('Number of positions = {0}\n'.format(len(positions)))
    if Optimizer.genealogy:
        Optimizer.Genealogyfile.write('\n')
    return

def update_outfile(ind, outfile):
    outfile.write('Individual {0}\n'.format(ind.index))
    outfile.write(repr(ind[0])+'\n')
    outfile.write('    Genealogy = {0}\n'.format(ind.history_index))
    outfile.write('    Energy = {0}\n'.format(ind.energy))
    outfile.write('    Fitness = {0}\n'.format(ind.fitness))
    outfile.write('    Swaplist = {0}\n'.format(ind.swaplist))

def update_structsumfile(ind, structsumfile):
    structsumfile.write(' Index = {0}\n'.format(ind.index))
    structsumfile.write('    Energy = {0}\n'.format(ind.energy))
    structsumfile.write('    Fitness = {0}\n'.format(ind.fitness))
    structsumfile.write('    Cell = {0}\n'.format(ind[0].get_cell()))
    structsumfile.write('    Pressure = {0}\n'.format(ind.pressure))
    structsumfile.write('    Genealogy = {0}\n'.format(ind.history_index))
    structsumfile.write('    Swaplist = {0}\n'.format(ind.swaplist))
    
def update_structfile(ind, structfile, Optimizer):
    if Optimizer.structure == 'Defect' or Optimizer.structure == 'Surface':
        sols = Atoms()
        sols.extend(ind[0])
        sols.extend(ind.bulki)
    elif Optimizer.structure == 'Crystal':
        sols = ind[0].repeat((3,3,3))
    else:
        sols = ind[0].copy()
    positions = sols.get_positions()
    if Optimizer.vacancy_output:
        for one in ind.vacancies:
            sols.append(Atom(symbol='X',position=one.position))
    Optimizer.output.write('Number of positions = {0}\n'.format(len(positions)))
    write_xyz(structfile, sols, ind.energy)
    return positions

def update_genealogy(ind, genefile):
    genefile.write('{0} '.format(ind.history_index))

def swaplist_check(ind, structype, outfile):
    sanchn = [[sym,0] for sym,c,m,u in Optimizer.atomlist]
    if structype=='Defect':
        solid = ind[0].copy()
        solid.extend(ind.bulki)
    else:
        solid = ind[0].copy()
    for i in range(len(sanchn)):
        nc = len([atm for atm in solid if atm.symbol==sanchn[i][0]])
        nc += [c for sym,c in ind.swaplist if sym==sanchn[i][0]][0]
        sanchn[i][1]=nc
        outfile.write('Atom structure {0} = {1}\n'.format(sanchn[i][0],
            sanchn[i][1]))