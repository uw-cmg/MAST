##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
from pymatgen.io.vaspio import *
import numpy as np
import logging
from MAST.utility.dirutil import *
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility import MASTObj
from MAST.utility import loggerutils
from pymatgen.core.structure import Structure
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.util.coord_utils import find_in_coord_list_pbc
from pymatgen.core.sites import PeriodicSite

class StructureExtensions(MASTObj):
    """Structure extensions
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'struc_work1': (Structure, None, 'First working Pymatgen Structure object (e.g. create a defect, or use work1 and work2 to interpolate positions)'),
            'struc_work2': (Structure, None, 'Second working Pymatgen Structure object'),
            'struc_init': (Structure, None, 'Initial structure at the beginning of the MAST recipe'),
            'name': (str, os.getenv("MAST_CONTROL"), 'Name of ingredient')
            }
        MASTObj.__init__(self, allowed_keys, **kwargs)
        self.logger = logging.getLogger(self.keywords['name'])
        self.logger = loggerutils.add_handler_for_recipe(self.keywords['name'], self.logger)

    def induce_defect(self, defect, coord_type, threshold):
        """Creates a defect, and returns the modified structure
            Args:
                defect <dict>: Defect subdictionary (single 
                               defect) of the form:
                        {'symbol': 'cr', 'type': 'interstitial', 
                    'coordinates': array([ 0. ,  0.5,  0. ])}}
                coord_type <str>: cartesian or fractional
                threshold <float>: Threshold for finding the
                                   defect position in what may
                                   be a relaxed, imperfect 
                                   structure.
            Returns:
                defected structure <Structure>
        """
        struct_ed = self.keywords['struc_work1'].copy() 
        symbol = defect['symbol'].title() #Cap first letter

        # If we have cartesian coordinates, then we convert them to fractional here.
        if ('cartesian' in coord_type):
            defect['coordinates'] = self._cart2frac(defect['coordinates'], self.keywords['struc_work1'])

        if (defect['type'] == 'vacancy'):
            self.logger.info('Creating a %s vacancy at %s' % (symbol, str(defect['coordinates'])))
            index = find_in_coord_list_pbc(self.keywords['struc_work1'].frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)
            if len(index) > 1:
                raise MASTError(self.__class__.__name__, "Multiple indices %s found. Check structure and/or adjust threshold %s to finer tolerance for ingredient %s" % (index, threshold,self.keywords['name']))
            if len(index) == 0:
                raise MASTError(self.__class__.__name__, "No indices found. Check structure and/or adjust threshold %s to lower tolerance for ingredient %s" % (threshold,self.keywords['name']))
            struct_ed.remove_sites([index])
        elif (defect['type'] == 'interstitial'):
            self.logger.info('Creating a %s interstitial at %s' % (symbol, str(defect['coordinates'])))
            struct_ed.append(symbol,
                                  defect['coordinates'],
                                  coords_are_cartesian=False,
                                  validate_proximity=True)
        elif (defect['type'] in ['antisite', 'substitution']):
            self.logger.info('Creating a %s antisite at %s' % (symbol, str(defect['coordinates'])))
            index = find_in_coord_list_pbc(self.keywords['struc_work1'].frac_coords,
                                       defect['coordinates'],
                                       atol=threshold)
            if len(index) > 1:
                raise MASTError(self.__class__.__name__, "Multiple indices %s found. Check structure and/or adjust threshold %s to finer tolerance for ingredient %s" % (index, threshold, self.keywords['name']))
            if len(index) == 0:
                raise MASTError(self.__class__.__name__, "No indices found. Check structure and/or adjust threshold %s to lower tolerance for ingredient %s" % (threshold, self.keywords['name']))
            struct_ed.replace(index, symbol)
            struct_ed = struct_ed.get_sorted_structure()
        else:
            raise RuntimeError('Defect type %s not supported' % defect['type'])

        return struct_ed
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
            1. Remove lines which closely match the neb moving 
               lines in the NEB section.
            2. Sort the structure by element and coordinate.
            3. Prepend the NEB moving lines to their element 
               sections.
            4. Return the modified structure.
            Args:
                folderstr <str>  : '00' = initial config, 
                                   '0N+1' = final config,
                                   '0N' = corresponding image
                neblines <list>  : list of NEB lines
                images <int>     : number of images
            Returns:
                sorted structure <Structure>
        """
        import MAST.data
        atol = 0.1 # Need fairly large tolerance to account for relaxation.
        sortedstruc = self.keywords['struc_work1'].get_sorted_structure()
        struct_ed = sortedstruc.copy()
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
                raise MASTError("pmgextend/structure_extensions", "No coordinate found matching %s for %s" % (mycoord, self.keywords['name']))
            nebidx.append(index[0]) #only take first site?
            mysite = sortedstruc.sites[index[0]]
            myelem = MAST.data.atomic_number[mysite.species_string]
            struct_ed.remove_sites([index[0]])
            struct_ed.insert(elemstarts[myelem], mysite.specie,
                                    mysite.frac_coords)
            sortedstruc = struct_ed.copy() #get new ordering
        if not len(nebidx) == len(neblines):
            raise MASTError("pmgextend/structure_extensions", "Not all NEB lines found for %s" % self.keywords['name'])
        return struct_ed

    def _get_element_indices(self, sortedstruc):
        """From a sorted structure, get the element indices
            Args:
                sortedstruc <Structure>: pymatgen structure 
                                sorted by electronegativity
            Returns:
                elstart <dict>: element index dictionary, in the
                                format:
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
        """Parse an NEB atomic movement line which has the 
           following format:
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
            Both struc_work1 and struc_work2 should be defined.
            Args:
                numim <int>: number of images
            Returns:
                <list of Structure>: list of interpolated 
                 structures, including the two endpoints
        """
        if self.keywords['struc_work1'] == None:
            raise MASTError(self.__class__.__name__, "No initial state structure for %s" % self.keywords['name'])
        if self.keywords['struc_work2'] == None:
            raise MASTError(self.__class__.__name__, "No final state structure for %s" % self.keywords['name'])
        structure_list = self.keywords['struc_work1'].interpolate(self.keywords['struc_work2'], numim+1)
        return structure_list

    def get_sd_array(self, phonon_center_site, phonon_center_radius, threshold=1e-1):
        """Create a selective dynamics array.
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius
                threshold <float>: absolute matching threshold for coordinates; 
                    e.g. if this value is 0.1 then both (0.1,0,0) and (0.3,0,0) 
                    will match for a search on center_site of (0.2,0,0)
        """
        mynbarr = self._get_neighbor_array(phonon_center_site, phonon_center_radius, self.keywords['struc_work1'], threshold)
        mysd = np.zeros([self.keywords['struc_work1'].num_sites,3],bool)
        for myn in mynbarr:
            mysd[myn]=np.ones(3,bool)
        return mysd

    def _get_neighbor_array(self, phonon_center_site, phonon_center_radius, mystruc, tol=1e-1):
        """
            Get a neighbor-index array.
            Use program_keywords 'phonon_center_site' and
            'phonon_center_radius' to limit the number of phonons calculated.
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius in Angstroms. If nonzero, all atoms in a radius around EACH site found in phonon_center_site will also be taken into account
                mystruc <Structure>: pymatgen Structure
                tol <float>: Tolerance for match-searching.
        """
        if phonon_center_site == None:
            return None
        pcscoord = np.array(phonon_center_site.strip().split(), float)
        tol = float(tol)
        pcsarr = find_in_coord_list_pbc(mystruc.frac_coords, pcscoord,tol)
        uniqsites = np.unique(pcsarr)

        if len(uniqsites) == 0:
            raise MASTError("pmgextend/structure_extensions", "No sites found for phonon centering for %s" % self.keywords['name'])

        if phonon_center_radius == None:
            return uniqsites

        nrad = float(phonon_center_radius)
        if nrad == 0:
            return uniqsites
        if nrad < 0:
            raise MASTError("pmgextend/structure_extensions", "Phonon center radius should not be less than zero for %s!" % self.keywords['name'])

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


    def get_multiple_sd_array(self, phonon_center_site, phonon_center_radius,threshold=1e-1):
        """Create a selective dynamics array, for use when every atom and every
            direction is a separate calculation (T F F, etc.)
            Args:
                phonon_center_site <str>: phonon center site (coordinate)
                phonon_center_radius <float>: phonon center radius, in Angstroms
                threshold <float>: absolute matching threshold for coordinates; e.g. if this value is 0.1, then both (0.1,0,0) and (0.3,0,0) will match for a search on center_site of (0.2,0,0)
            Returns:
                mysdlist <list>: list of SD arrays
        """
        if phonon_center_site == None:
            return None
        mynbarr = self._get_neighbor_array(phonon_center_site, phonon_center_radius, self.keywords['struc_work1'], threshold)
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
                coordstruc <Structure>: Structure object with
                    the coordinates for grafting
                self.keywords['struc_work1'] will contain
                    the elements and lattice parameter, which
                    will not be touched.
            Returns:
                modified Structure object <Structure>
        """
        goodstruc = self.keywords['struc_work1'].copy()
        lengoodsites=len(goodstruc.sites)
        lencoordsites=len(coordstruc.sites)
        if not (lengoodsites == lencoordsites):
            raise MASTError(self.__class__.__name__, "Original and coordinate structures do not have the same amount of sites in %s" % self.keywords['name'])
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
    def strain_lattice(self, strainstring):
        """Strain the lattice. 
            Args:
                strainstring <str>: Strain string from 
                    'mast_strain' in program keys, 
                    with space-separated strain along 
                    lattice vectors a, b, and c, e.g.
                    "0.98 0.96 1.01"
            Returns:
                newstructure <Structure>: strained structure
        """
        mystructure = self.keywords['struc_work1']
        strainsplit = strainstring.strip().split()
        strarray = np.array([[0.],[0.],[0.]],'float')
        for sidx in range(0,3):
            strflt = float(strainsplit[sidx])
            strarray[sidx][0]=strflt
        newlattice = Lattice(np.multiply(mystructure._lattice.matrix, strarray)) #be very careful here. np.multiply is NOT regular matrix multiplication.
        newstructure = mystructure.copy()
        newstructure.modify_lattice(newlattice)
        return newstructure

    def scale_structure(self, scale):
        """Scale the structure.
            Args:
                oldstr <structure>: old structure to be scaled
                scale <str>: scale parameter; should be >= 2
                        Should be either a single digit "2"
                        or a string of 3 digits "1 2 3"
                        signifying scaling by 1a x 2b x 3c
            Returns:
                scaledstr <structure>
        """
        scaledstr = self.keywords['struc_work1'].copy()
        if type(scale) == int:
            scaleinput = scale
        else:
            scale = scale.strip()
            scalesplit = scale.split()
            if len(scalesplit) == 1:
                scaleinput = float(scalesplit[0])
            elif len(scalesplit) == 3:
                scaleinput = map(float, scalesplit)
            else:
                self.logger.error("Wrong number of inputs for scaling: %s " % scalesplit)
                raise MASTError("Wrong number of inputs for scaling: %s " % scalesplit)
                return None
        scaledstr.make_supercell(scaleinput)
        return scaledstr

    def scale_defect(self, defect, coord_type, threshold):
        """Scales the defect dictionary and returns the modified structure
            Args:
                keyword struc_work1: scaled-up structure
                keyword struc_work2: original structure
                defect <dict>: Defect subdictionary (single 
                               defect) of the form:
                        {'symbol': 'cr', 'type': 'interstitial', 
                    'coordinates': array([ 0. ,  0.5,  0. ])}}
                coord_type <str>: cartesian or fractional
                threshold <float>: Threshold for finding the
                                   defect position in what may
                                   be a relaxed, imperfect 
                                   structure.
            Returns:
                defected structure <Structure>
        """
        mycoords = defect['coordinates']
        origstr = self.keywords['struc_work2']
        scaledstr = self.keywords['struc_work1']
        olda = origstr.lattice.a
        oldb = origstr.lattice.b
        oldc = origstr.lattice.c
        newa = scaledstr.lattice.a
        newb = scaledstr.lattice.b
        newc = scaledstr.lattice.c
        #if newa is twice as big as olda, coordinate should be half as big
        coorda = mycoords[0]*(olda/newa)
        coordb = mycoords[1]*(oldb/newb)
        coordc = mycoords[2]*(oldc/newc)
        newdict = dict(defect)
        newdict['coordinates'] = np.array([coorda, coordb, coordc],'float')
        returnstr = self.induce_defect(newdict, coord_type, threshold)
        return returnstr
    def scale_defect_by_LMN(self, scale, defect, coord_type, threshold):
        """Scales the defect dictionary and returns the modified structure
            Args:
                keyword struc_work1: scaled-up structure
                scale <str>: String "LxMxN" for scaling defect position
                defect <dict>: Defect subdictionary (single 
                               defect) of the form:
                        {'symbol': 'cr', 'type': 'interstitial', 
                    'coordinates': array([ 0. ,  0.5,  0. ])}}
                coord_type <str>: cartesian or fractional
                threshold <float>: Threshold for finding the
                                   defect position in what may
                                   be a relaxed, imperfect 
                                   structure.
            Returns:
                defected structure <Structure>
        """
        mycoords = defect['coordinates']
        scaledstr = self.keywords['struc_work1']
        scalelist = scale.strip().split("x")
        scale1 = float(scalelist[0])
        scale2 = float(scalelist[1])
        scale3 = float(scalelist[2])
        #if newa is twice as big as olda, coordinate should be half as big
        coorda = mycoords[0]/scale1
        coordb = mycoords[1]/scale2
        coordc = mycoords[2]/scale3
        newdict = dict(defect)
        newdict['coordinates'] = np.array([coorda, coordb, coordc],'float')
        returnstr = self.induce_defect(newdict, coord_type, threshold)
        return returnstr
