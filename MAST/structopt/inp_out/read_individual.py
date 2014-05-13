from MAST.structopt.inp_out.read_xyz import read_xyz
from MAST.structopt.generate.Individual import Individual
from ase import Atom, Atoms

def read_individual(indivfile, n=-1):
    """Function to write the data of an individual class object to a flat file
    Input:
        indivfile = String or fileobject for file to be read from
        n = which individual from file to return. Default is last individual written.
            optional All
    Output:
        returns an individual class object or list of individual class objects depending on value of n
    """
    if isinstance(indivfile, str):
        indivfile=open(indivfile, 'r')
    all_lines = indivfile.readlines()
    indivfile.close()
    linen = 0
    all_indivs = []
    while linen < len(all_lines):
        if '----------' in all_lines[linen]:
            individ = Individual(Atoms())
        elif 'Structure information' in all_lines[linen]:
            natomstruct = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natomstruct):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ = Individual(atomstruct)
            linen += 2+natomstruct
        elif 'structure cell' in all_lines[linen]:
            cell_line = all_lines[linen].split('=')
            structcell = eval(cell_line[1])
            individ[0].set_cell(structcell)
        elif 'fitness' in all_lines[linen]:
            fitline = all_lines[linen].split('=')
            individ.fitness = float(fitline[1])
        elif 'history_index' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.history_index = line[1].strip()
        elif 'index' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.index = float(line[1])
        elif 'tenergymx' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.tenergymx = float(line[1])
        elif 'tenergymin' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.tenergymin = float(line[1])
        elif 'energy' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.energy = float(line[1])
        elif 'pressure' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.pressure = float(line[1])
        elif 'volume' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.volume = float(line[1])
        elif 'force' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.force = float(line[1])
        elif 'purebulkenpa' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.purebulkenpa = float(line[1])
        elif 'natomsbulk' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.natomsbulk = float(line[1])
        elif 'fingerprint' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.fingerprint = eval(line[1])
        elif 'swpalist' in all_lines[linen]:
            line = all_lines[linen].split('=')
            individ.swaplist = eval(line[1])
        elif 'bulki cell' in all_lines[linen]:
            line = all_lines[linen].split('=')
            cell = eval(line[1])
            individ.bulki.set_cell(cell)
        elif 'bulki' in all_lines[linen]:
            natoms = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natoms):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ.bulki = atomstruct.copy()
            linen += 2+natoms
        elif 'bulko cell' in all_lines[linen]:
            line = all_lines[linen].split('=')
            cell = eval(line[1])
            individ.bulko.set_cell(cell)
        elif 'bulko' in all_lines[linen]:
            natoms = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natoms):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ.bulko = atomstruct.copy()
            linen += 2+natoms
        elif 'box cell' in all_lines[linen]:
            line = all_lines[linen].split('=')
            cell = eval(line[1])
            individ.box.set_cell(cell)
        elif 'box' in all_lines[linen]:
            natoms = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natoms):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ.box = atomstruct.copy()
            linen += 2+natoms
        elif 'vacancies cell' in all_lines[linen]:
            line = all_lines[linen].split('=')
            cell = eval(line[1])
            individ.vacancies.set_cell(cell)
        elif 'vacancies' in all_lines[linen]:
            natoms = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natoms):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ.vacancies = atomstruct.copy()
            linen += 2+natoms
        elif 'swaps cell' in all_lines[linen]:
            line = all_lines[linen].split('=')
            cell = eval(line[1])
            individ.swaps.set_cell(cell)
        elif 'swaps' in all_lines[linen]:
            natoms = int(all_lines[linen+1])
            atomstruct = Atoms()
            for i in range(natoms):
                a = all_lines[linen+i+3].split()
                sym = a[0]
                position = [float(a[1]),float(a[2]),float(a[3])]
                atomstruct.append(Atom(symbol=sym,position=position))
            individ.swaps = atomstruct.copy()
            linen += 2+natoms
        elif 'Finish' in all_lines[linen]:
            all_indivs.append(individ.duplicate())
        linen+=1
    if n=='All':
        return all_indivs
    else:
        return all_indivs[n]
    
    