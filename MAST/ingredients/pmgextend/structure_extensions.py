from pymatgen.io.vaspio import *
import numpy as np
from MAST.utility.dirutil import *
from MAST.utility import MASTError
from MAST.utility import MASTFile

def induce_defect(base_structure, defect, coord_type, threshold):
    """Creates a defect, and returns the modified structure
        mast_defect is a dictionary like this: 
        'defect_1': {'symbol': 'cr', 'type': 'interstitial', 
                    'coordinates': array([ 0. ,  0.5,  0. ])}}
        'defect_2': {'symbol': 'o', 'type': 'vacancy', 
                    'coordinates': array([ 0.25,  0.25,  0.25])}
        'coord_type': 'fractional' 
    """
    #print 'Defect in induce_defect', defect
    #print 'base_structure in induce_defect', base_structure
    struct_ed = StructureEditor(base_structure) #should be updated using get_new_structure)
    symbol = defect['symbol'].title() #Cap first letter

    # If we have cartesian coordinates, then we convert them to fractional here.
    if ('cartesian' in coord_type):
        defect['coordinates'] = _cart2frac(defect['coordinates'],base_structure)

    if (defect['type'] == 'vacancy'):
        print 'Creating a %s vacancy at %s' % (symbol, str(defect['coordinates']))

        #print defect['coordinates']
        index = find_in_coord_list(base_structure.frac_coords,
                                   defect['coordinates'],
                                   atol=threshold)
        #print base_structure.frac_coords
        #print 'Index of deleted atom is', index
        struct_ed.delete_site(index)
    elif (defect['type'] == 'interstitial'):
        print 'Creating a %s interstitial at %s' % (symbol, str(defect['coordinates']))

        struct_ed.append_site(symbol,
                              defect['coordinates'],
                              coords_are_cartesian=False,
                              validate_proximity=True)
    elif (defect['type'] in ['antisite', 'substitution']):
        print 'Creating a %s antisite at %s' % (symbol, str(defect['coordinates']))

        index = find_in_coord_list(base_structure.frac_coords,
                                   defect['coordinates'],
                                   atol=threshold)

        struct_ed.replace_site(index, symbol)
    else:
        raise RuntimeError('Defect type %s not supported' % defect['type'])

    return struct_ed.modified_structure
def _cart2frac(position, base_structure):
    """Converts between cartesian coordinates and fractional coordinates"""
    fractional =  base_structure.lattice.get_fractional_coords(position)
    for i in range(len(fractional)):
        if (fractional[i] < 0.0):
            fractional[i] += 1.0
        elif (fractional[i] > 1.0):
            fractional[i] -= 1.0
    return fractional

