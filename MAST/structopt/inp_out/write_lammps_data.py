import os
from MAST.structopt.tools.lammps import prism

def write_lammps_data(data_name, atoms):
    """Function to write ase atoms object to LAMMPS file
    Input:
        data_name = String or fileobject for data to be written to
        atoms = ase atoms object to be written
    Output:
        None. Data file will be written
    """
    if isinstance(data_name, str):
        f = open(data_name, 'w')
    else:
        f = data_name
    f.write('\n\n')
    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    f.write('{0} \t atoms \n'.format(n_atoms))
    #Order species alphabetically
    ordered_symbols = sorted(list(set(symbols)))
    #Write these symbols to a file to ensure they are kept
    dir = os.path.dirname(data_name)
    afilename = os.path.join(dir, 'atom_symbols')
    a = open(afilename,'w')
    for sym in ordered_symbols:
        a.write('{0}\n'.format(sym))
    a.close()
    n_types = len(ordered_symbols)
    f.write('{0}  atom types\n'.format(n_types))
    p = prism(atoms.get_cell())
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()

    f.write('0.0 %s  xlo xhi\n' % xhi)
    f.write('0.0 %s  ylo yhi\n' % yhi)
    f.write('0.0 %s  zlo zhi\n' % zhi)

    if p.is_skewed():
        f.write('%s %s %s  xy xz yz\n' % (xy, xz, yz))
    f.write('\n\n')

    f.write('Atoms \n\n')
    for i, r in enumerate(map(p.pos_to_lammps_str,
                              atoms.get_positions())):
        s = ordered_symbols.index(symbols[i]) + 1
        f.write('%6d %3d %s %s %s\n' % ((i+1, s)+tuple(r)))
    f.close()

