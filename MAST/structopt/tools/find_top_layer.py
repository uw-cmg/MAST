from ase import Atom, Atoms

def find_top_layer(surf,topthick):
    """Development function for identifying viable surface atoms.
    *** needs development ***
    """
    top=Atoms()
    bulk=Atoms()
    zs=[z for x,y,z in surf.get_positions()]
    thick=max(zs)-min(zs)
    for one in surf:
        x,y,z=one.position
        if z>=max(zs)-topthick:
            top.append(one)
        else:
            bulk.append(one)
    return top, bulk

