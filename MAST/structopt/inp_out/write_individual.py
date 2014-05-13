from MAST.structopt.inp_out.write_xyz import write_xyz

def write_individual(individ, indivfile):
    """Function to write the data of an individual class object to a flat file
    Input:
        individ = Individual class object to be written
        indivfile = String or fileobject for file to be written to 
    Output:
        No output returned.  Information is written to file
    """
    if isinstance(indivfile, str):
        indivfile=open(indivfile, 'a')
    #Write break
    indivfile.write('----------\n')
    #Write structure information
    indivfile.write('Structure information\n')
    write_xyz(indivfile, individ[0])
    indivfile.write('structure cell = {0}\n'.format(get_atom_cell(individ[0])))
    #Write additional information
    indivfile.write('fitness = {0}\n'.format(individ.fitness))
    indivfile.write('index = {0}\n'.format(individ.index))
    indivfile.write('history_index = {0}\n'.format(individ.index))
    indivfile.write('energy = {0}\n'.format(individ.energy))
    indivfile.write('tenergymx = {0}\n'.format(individ.tenergymx))
    indivfile.write('tenergymin = {0}\n'.format(individ.tenergymin))
    indivfile.write('pressure = {0}\n'.format(individ.pressure))
    indivfile.write('volume = {0}\n'.format(individ.volume))
    indivfile.write('force = {0}\n'.format(individ.force))
    indivfile.write('purebulkenpa = {0}\n'.format(individ.purebulkenpa))
    indivfile.write('natomsbulk = {0}\n'.format(individ.natomsbulk))
    indivfile.write('fingerprint = {0}\n'.format(individ.fingerprint))
    indivfile.write('swaplist = {0}\n'.format(individ.swaplist))
    #Write additional structure information
    indivfile.write('bulki\n')
    write_xyz(indivfile, individ.bulki)
    indivfile.write('bulki cell = {0}\n'.format(get_atom_cell(individ.bulki)))
    indivfile.write('bulko\n')
    write_xyz(indivfile, individ.bulko)
    indivfile.write('bulko cell = {0}\n'.format(get_atom_cell(individ.bulko)))
    indivfile.write('box\n')
    write_xyz(indivfile, individ.box)
    indivfile.write('box cell = {0}\n'.format(get_atom_cell(individ.box)))
    indivfile.write('vacancies\n')
    write_xyz(indivfile, individ.vacancies)
    indivfile.write('vacancies cell = {0}\n'.format(get_atom_cell(individ.vacancies)))
    indivfile.write('swaps\n')
    write_xyz(indivfile, individ.swaps)
    indivfile.write('swaps cell = {0}\n'.format(get_atom_cell(individ.swaps)))
    indivfile.write('Finish')
    indivfile.close()
    return

def get_atom_cell(atomsobj):
    cell = atomsobj.get_cell()
    cell_list = []
    for i in range(3):
        clist = []
        for j in range(3):
            clist.append(cell[i][j])
        cell_list.append(clist)
    return cell_list