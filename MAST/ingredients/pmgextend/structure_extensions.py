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

def sort_structure_and_neb_lines(mystruc, basestruc, neblines, initorfin):
    """Sort the structure in an approved way:
        2. Remove lines which closely match the neb moving lines in the
        NEB section.
        3. Sort the structure by element and coordinate.
        4. Prepend the NEB moving lines to their element sections.
        5. Return the modified structure.

        mystruc <Structure>   : pymatgen structure of endpoint
        basestruc <Structure> : pymatgen structure of initial 
                                cell (used to match vacancies)
        initorfin <int>       : 0 = initial cfg, 1 = final cfg
        neblines <list>       : list of NEB lines
    """
    import MAST.data
    atol = 0.1 # Need fairly large tolerance to account for relaxation.
    sortedstruc = mystruc.get_sorted_structure()
    sortedstruc_base = basestruc.get_sorted_structure()
    struct_ed = StructureEditor(sortedstruc)
    struct_ed_base = StructureEditor(sortedstruc_base)
    # We are not actually translating the sites; just get them into unit
    # cell
    #TTM DEBUG REMOVE THIS struct_ed.translate_sites(range(0,len(sortedstruc.sites)),np.zeros(3),True)
    nebidx = list()

    elemstarts = get_element_indices(sortedstruc)
    for nebline in neblines:
        usebase=0
        nebdict = _parse_neb_line(nebline)
        mycoord = nebdict['coord'][initorfin]
        index = find_in_coord_list(sortedstruc.frac_coords, mycoord, atol)
        print 'TTM DEBUG: index: ', index
        if len(index) == 0: #try the base structure, for vacancies
            usebase=1
            index = find_in_coord_list(sortedstruc_base.frac_coords,
                mycoord, atol)
        if len(index) == 0:
            raise MASTError("pmgextend/structure_extensions", "No coordinate found matching %s" % mycoord)
        if usebase==0:
            nebidx.append(index[0]) #only take first site?
            mysite = sortedstruc.sites[index[0]]
            myelem = MAST.data.atomic_number[mysite.species_string]
            struct_ed.delete_site(index)
            struct_ed.insert_site(elemstarts[myelem], mysite.specie,
                                    mysite.frac_coords)
        else:
            nebidx.append(index[0]) #only take first site?
            mysite = sortedstruc_base.sites[index[0]]
            myelem = MAST.data.atomic_number[mysite.species_string]
            struct_ed_base.delete_site(index)
            struct_ed_base.insert_site(elemstarts[myelem], mysite.specie,
                                    mysite.frac_coords)

    if not len(nebidx) == len(neblines):
        raise MASTError("pmgextend/structure_extensions", "Not all NEB lines found.")
    return struct_ed.modified_structure

def get_element_indices(sortedstruc):
    """From a sorted structure, get the element indices
        Args:
            sortedstruc <Structure>: pymatgen structure sorted by
                                        electronegativity
        Returns:
            elstart <dict>: element index dictionary, in the format:
                        elstart[element number] = starting index, e.g.
                        elstart[24] = 0
    """
    elstart = dict()
    numlist = sortedstruc.atomic_numbers
    listlen = len(numlist)
    idx = 0
    while idx < listlen:
        atomnum = numlist[idx]
        if not atomnum in elstart.keys():
            elstart[atomnum] = idx
        idx = idx + 1
    return elstart

def _parse_neb_line(nebline):
    """Parse an NEB atomic movement line which has the following format:
                Cr, 0 0 0, 0.5 0.5 0.5 in the input file, and
                ['Cr','0 0 0','0.5 0.5 0.5'] here as nebline
                element, init_coord, fin_coord
        Args:
            line <list of str>: NEB atomic movement line
        Returns:
            linedict <dict>:
                             ['element'] <int> = element's atomic number
                             ['coord'][0] <np.array> = initial coordinates
                             ['coord'][1] <np.array> = final coordinates
    """
    nebdict=dict()
    #print "TTM DEBUG: nebline", nebline
    nebelem = str(nebline[0].strip()).title()
    import MAST.data
    nebdict['element'] = MAST.data.atomic_number[nebelem]
    nebdict['coord'] = dict()
    nebdict['coord'][0] = np.array(nebline[1].split(), dtype='float')
    nebdict['coord'][1] = np.array(nebline[2].split(), dtype='float')
    return nebdict

def do_interpolation(parentstructures, numim):
    """Do interpolation.
        Args:
            parentstructures <list of Structure>: list of endpoint
                        pymatgen Structure objects
            numim <int>: number of images
    """
    if parentstructures == None:
        raise MASTError("pmgextend/structure_extensions","Bad number of parent paths.")
    struct_init = parentstructures[0]
    struct_fin = parentstructures[1]
    structure_list = struct_init.interpolate(struct_fin, numim+1) + 1 )
    return structure_list

