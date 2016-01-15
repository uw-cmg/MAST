import numpy

def shift_atoms(atoms, position=None):
    """Function to shift the positions of the atoms so that
    the given location is in the center of the cell.
    Input:
        atoms = ASE atoms object to be shifted
        position = [x,y,z] position to shift to center of cell
            default is center of mass of atoms
    Output:
        atoms = ASE atoms object that has been shifted
    """
    atoms = atoms.copy()
    if not position:
        position = atoms.get_center_of_mass()
    cell = numpy.maximum.reduce(atoms.get_cell())
    trans = [cell[i]/2.0-position[i] for i in range(3)]
    atoms.translate(trans)
    positions = atoms.get_positions()
    for i in range(len(positions)):
        if atoms.pbc[0]:
            if trans[0] > 0:
                if positions[i][0] > cell[0]:
                    positions[i][0] -= cell[0]
            else:
                if positions[i][0] < 0:
                    positions[i][0] += cell[0]
        if atoms.pbc[1]:
            if trans[1] > 0:
                if positions[i][1] > cell[1]:
                    positions[i][1] -= cell[1]
            else:
                if positions[i][1] < 0:
                    positions[i][1] += cell[1]
        if atoms.pbc[2]:
            if trans[2] >0:
                if positions[i][2] > cell[2]:
                    positions[i][2] -= cell[2]
            else:
                if positions[i][2] < 0:
                    positions[i][2] += cell[2]
    atoms.set_positions(positions)
    return atoms

