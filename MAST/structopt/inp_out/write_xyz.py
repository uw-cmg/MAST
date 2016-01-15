from ase import Atom, Atoms

def write_xyz(fileobj,atms,data=0):
    """Function to write xyz file with some data
    adapted from ase.io.xyz"""
    if isinstance(fileobj, str):
        fileobj=open(fileobj, 'a')
    symbols = atms.get_chemical_symbols()
    natoms = len(symbols)
    fileobj.write('%d\n' % natoms)
    fileobj.write(repr(data)+'\n')
    for s, (x, y, z) in zip(symbols, atms.get_positions()):
        fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
    return
