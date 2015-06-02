try:
    from ase import Atom, Atoms
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt write_xyz.py, ASE must be installed."

def write_xyz(fileobj,atms,data=0):
    """Function to write xyz file with some data
    adapted from ase.io.xyz"""
    if isinstance(fileobj, str):
        fileobj=open(fileobj, 'a')
    symbols = atms.get_chemical_symbols()
    natoms = len(symbols)
    fileobj.write('%d\n' % natoms)
    #fileobj.write(repr(data)+'\n')
    cell = atms.get_cell()
    fileobj.write('%s %10.5f %10.5f %10.5f \n' % (repr(data),cell[0][0],cell[1][1],cell[2][2]))
    for s, (x, y, z) in zip(symbols, atms.get_positions()):
        fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
    return
