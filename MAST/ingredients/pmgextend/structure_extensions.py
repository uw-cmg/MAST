from pymatgen.io.vaspio import *
import numpy as np
from MAST.utility.dirutil import *
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility import MASTObj
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.util.coord_utils import find_in_coord_list_pbc
from pymatgen.core.sites import PeriodicSite
import logging

class StructureExtensions(MASTObj):
    """Structure extensions
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'struc_work1': (Structure, None, 'First working Pymatgen Structure object (e.g. create a defect, or use work1 and work2 to interpolate positions)'),
            'struc_work2': (Structure, None, 'Second working Pymatgen Structure object'),
            'struc_init': (Structure, None, 'Initial structure at the beginning of the MAST recipe')
            }
        MASTObj.__init__(self, allowed_keys, **kwargs)

        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)

    def induce_defect(self, defect, coord_type, threshold):
        """Creates a defect, and returns the modified structure
            mast_defect is a dictionary like this: 
            'defect_1': {'symbol': 'cr', 'type': 'interstitial', 
                        'coordinates': array([ 0. ,  0.5,  0. ])}}
            'defect_2': {'symbol': 'o', 'type': 'vacancy', 
                        'coordinates': array([ 0.25,  0.25,  0.25])}
            'coord_type': 'fractional' 
        """
        struct_ed = StructureEditor(self.keywords['struc_work1']) #should be updated using get_new_structure)
        symbol = defect['symbol'].title() #Cap first letter

        # If we have cartesian coordinates, then we convert them to fractional here.
        if ('cartesian' in coord_type):
            defect['coordinates'] = self._cart2frac(defect['coordinates'], self.keywords['struc_work1'])

        if (defect['type'] == 'vacancy'):
            self.logger.info('Creating a %s vacancy at %s' % (symbol, str(defect['coordinates'])))

            index = find_in_coord_list(self.keywords['struc_work1'].frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)
            struct_ed.delete_site(index)
        elif (defect['type'] == 'interstitial'):
            self.logger.info('Creating a %s interstitial at %s' % (symbol, str(defect['coordinates'])))

            struct_ed.append_site(symbol,
                                  defect['coordinates'],
                                  coords_are_cartesian=False,
                                  validate_proximity=True)
        elif (defect['type'] in ['antisite', 'substitution']):
            self.logger.info('Creating a %s antisite at %s' % (symbol, str(defect['coordinates'])))

            index = find_in_coord_list(self.keywords['struc_work1'].frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)

            struct_ed.replace_site(index, symbol)
        else:
            raise RuntimeError('Defect type %s not supported' % defect['type'])

        return struct_ed.modified_structure
    def _cart2frac(self, position, base_structure):
        """Converts between cartesian coordinates and fractional coordinates"""
        fractional =  base_structure.lattice.get_fractional_coords(position)
        for i in range(len(fractional)):
            if (fractional[i] < 0.0):
                fractional[i] += 1.0
            elif (fractional[i] > 1.0):
                fractional[i] -= 1.0
        return fractional

    def sort_structure_and_neb_lines(self, neblines, folderstr, images=0):
        """Sort the structure in an approved way:
            2. Remove lines which closely match the neb moving lines in the
            NEB section.
            3. Sort the structure by element and coordinate.
            4. Prepend the NEB moving lines to their element sections.
            5. Return the modified structure.

            folderstr <str>       : '00' = initial config, '0N+1' = final config,
                                    '0N' = corresponding image
            neblines <list>       : list of NEB lines
            images <int>          : number of images
        """
        import MAST.data
        atol = 0.1 # Need fairly large tolerance to account for relaxation.
        sortedstruc = self.keywords['struc_work1'].get_sorted_structure()
        struct_ed = StructureEditor(sortedstruc)
        nebidx = list()
        elemstarts = self._get_element_indices(sortedstruc)
        for nebline in neblines:
            nebdict = self._parse_neb_line(nebline)
            temp_fin = sortedstruc.copy()
            temp_fin.append(nebdict['element'],nebdict['coord'][1])
            temp_start = sortedstruc.copy()
            temp_start.append(nebdict['element'],nebdict['coord'][0])
            strlist=temp_start.interpolate(temp_fin, images+1)
            lastidx = strlist[0].num_sites-1
            if folderstr == '00':
                mycoord = nebdict['coord'][0]
            elif folderstr == str(images+1).zfill(2):
                mycoord = nebdict['coord'][1]
            else:
                mystridx = int(folderstr)
                mycoord = strlist[mystridx].frac_coords[lastidx]
            
            indexraw = find_in_coord_list_pbc(sortedstruc.frac_coords, mycoord, atol)
            index=list()
            for indexentry in indexraw:
                if sortedstruc.species[indexentry] == temp_start.species[lastidx]:
                    index.append(indexentry)
            if len(index) == 0:
                raise MASTError("pmgextend/structure_extensions", "No coordinate found matching %s" % mycoord)
            nebidx.append(index[0]) #only take first site?
            mysite = sortedstruc.sites[index[0]]
            myelem = MAST.data.atomic_number[mysite.species_string]
            struct_ed.delete_site(index[0])
            struct_ed.insert_site(elemstarts[myelem], mysite.specie,
                                    mycoord)
        if not len(nebidx) == len(neblines):
            raise MASTError("pmgextend/structure_extensions", "Not all NEB lines found.")
        return struct_ed.modified_structure

    def _get_element_indices(self, sortedstruc):
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

    def _parse_neb_line(self, nebline):
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
        nebelem = str(nebline[0].strip()).title()
        import MAST.data
        nebdict['element'] = MAST.data.atomic_number[nebelem]
        nebdict['coord'] = dict()
        nebdict['coord'][0] = np.array(nebline[1].split(), dtype='float')
        nebdict['coord'][1] = np.array(nebline[2].split(), dtype='float')
        return nebdict


    def do_interpolation(self, numim):
        """Do interpolation.
            Args:
                numim <int>: number of images
        """
        if self.keywords['struc_work1'] == None:
            raise MASTError(self.__class__.__name__, "No initial state structure")
        if self.keywords['struc_work2'] == None:
            raise MASTError(self.__class__.__name__, "No final state structure")
        structure_list = self.keywords['struc_work1'].interpolate(self.keywords['struc_work2'], numim+1)
        return structure_list

    def get_sd_array(self, phonon_center_site, phonon_center_radius):
        """Create a selective dynamics array.
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius
        """
        mynbarr = self._get_neighbor_array(phonon_center_site, phonon_center_radius, self.keywords['struc_work1'])
        mysd = np.zeros([self.keywords['struc_work1'].num_sites,3],bool)
        for myn in mynbarr:
            mysd[myn]=np.ones(3,bool)
        return mysd

    def _get_neighbor_array(self, phonon_center_site, phonon_center_radius, mystruc, tol=1e-1):
        """
            Get a neighbor-index array.
            Use program_keywords 'phonon_center_site' and
            'phonon_center_radius' to limit the number of phonons calculated.
            ['program_keys']['phonon'][label]['phonon_center_site']
                    should be a coordinate
                    If the key is missing, all atoms will be taken into account.
            ['program_keys']['phonon'][label]['phonon_center_radius']
                    should be a positive float in ANGSTROMS (Not fractional.)
                    If the key is missing or 0, nothing extra happens.
                    If the key is present and nonzero, then all atoms in a
                        radius around EACH site found in phonon_center_site
                        will also be taken into account.
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius
                mystruc <Structure>: pymatgen Structure
                tol <float>: Tolerance for match-searching.
        """
        if phonon_center_site == None:
            return None
        pcscoord = np.array(phonon_center_site.strip().split(), float)
        pcsarr = find_in_coord_list(mystruc.frac_coords, pcscoord,tol)
        uniqsites = np.unique(pcsarr)

        if len(uniqsites) == 0:
            raise MASTError("pmgextend/structure_extensions", "No sites found for phonon centering.")

        if phonon_center_radius == None:
            return uniqsites

        nrad = float(phonon_center_radius)
        if nrad == 0:
            return uniqsites
        if nrad < 0:
            raise MASTError("pmgextend/structure_extensions", "Phonon center radius should not be less than zero!")

        nbtotarr=None
        for pcs in uniqsites:
            neighbors = mystruc.get_neighbors(mystruc[pcs], nrad, True)
            if nbtotarr == None:
                nbtotarr = neighbors
            else:
                np.concatenate([nbtotarr, neighbors])
        nbsitelist=list()
        for nbr in nbtotarr:
            nbsitelist.append(nbr[-1])
        nbsitelist = np.array(nbsitelist)
        alltotarr = np.concatenate([uniqsites, nbsitelist])
        allsites = np.unique(alltotarr)
        return allsites


    def get_multiple_sd_array(self, phonon_center_site, phonon_center_radius):
        """Create a selective dynamics array, for use when every atom and every
            direction is a separate calculation (T F F, etc.)
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius, in Angstroms
            Returns:
                mysdlist <list>: list of SD arrays
        """
        if phonon_center_site == None:
            return None
        mynbarr = self._get_neighbor_array(phonon_center_site, phonon_center_radius, self.keywords['struc_work1'])
        mysdlist=list()
        for myn in mynbarr:
            for myct in range(0,3):
                mysd = np.zeros([self.keywords['struc_work1'].num_sites,3],bool)
                mysd[myn]=np.zeros(3,bool)
                mysd[myn][myct]=1
                mysdlist.append(mysd)
        return mysdlist

    def graft_coordinates_onto_structure(self, coordstruc):
        """Graft coordinates from mast_coordinates Structure objects
            onto the appropriate structure
            Args:
                coordstrucs <list>: Structure object with
                    the coordinates for grafting
                self.keywords['struc_work1'] will contain
                    the elements and lattice parameter, which
                    will not be touched.
            Returns:
                modstruc: modified Structure objects
        """
        goodstruc = self.keywords['struc_work1'].copy()
        lengoodsites=len(goodstruc.sites)
        lencoordsites=len(coordstruc.sites)
        if not (lengoodsites == lencoordsites):
            raise MASTError(self.__class__.__name__, "Original and coordinate structures do not have the same amount of sites.")
        cct=0
        newsites=list()
        mylattice=goodstruc.lattice
        while cct < lengoodsites:
            newcoords=coordstruc.sites[cct].frac_coords
            oldspecie=goodstruc.sites[cct].specie
            newsite=PeriodicSite(oldspecie, newcoords, mylattice)
            newsites.append(newsite)
            cct=cct+1
        goodstruc.remove_sites(range(0,lengoodsites))
        for cct in range(0, lengoodsites):
            goodstruc.append(newsites[cct].specie,
                newsites[cct].frac_coords)
        return goodstruc

