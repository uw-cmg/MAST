def check_cell_type(indiv):
    """Basic check of cell type
    *** needs development***
    """
    cell=indiv.get_cell()
    a,b,c=cell
    if b[0] != 0:
        if c[0] != 0:
            cell_type='triclinic'
        else:
            cell_type='hexagonal'
    elif c[0] != 0:
        cell_type='monoclinic'
    elif b[1] != c[2]:
        if b[1] != a[0]:
            cell_type='orthorhombic'
        else:
            cell_type='tetragonal'
    else:
        cell_type='cubic'
    return cell_type

