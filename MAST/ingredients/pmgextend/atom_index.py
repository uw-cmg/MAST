##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2015-02-12 by Tam Mayeshiba
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
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions as SE
from MAST.utility import InputOptions

class AtomIndex(MASTObj):
    """Atom Index
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            #'startSE': (SE, None, 'Structure extensions class')
            'input_options': (InputOptions, None, 'Input options')
    }
    self.startSE = ""
    self.startdict = ""
    self.startdefects = ""
    self.startnebs = ""
    self.scalingSEs = dict()
    self.scalingdicts = dict()
    self.scalingdefects = dict()
    self.scalingnebs = ""
    self.input_options = self.keywords['input_options']
    return

    def set_up_initial_index(self):
        """Set up the initial index (first time).
        """
        self.startSE = self.get_start_structure()
        self.startdict = self.build_structure_dictionary(self, self.startSE)


    def get_start_structure(self):
        """Get starting structure. 
        """
        startstr=self.input_options.get_item('structure','structure')
        return SE(struc_work1 = startstr)
    

    def build_structure_dictionary(self, mySE):
        """Build a structure dictionary.
        """
        sdict=dict()
        sdx=0
        for site in mySE.keywords['struc_work1'].sites:
            skey=self.convert_int_to_atomidx(sdx)
            sdict[skey]=dict()
            sdict[skey]['original_frac_coords']=site.frac_coords
            sdict[skey]['element']=site.species_string
            sdict[skey]['specie']=site.specie
            sdx=sdx+1
        return sdict
    
    def make_entry_from_site(self):
        """
        """
        return

    def add_atom_specific_keywords_to_structure_dictionary(self):
        """Add atom specific keywords using the coordinates section
        """
        return
    
    def add_element_specific_keywords_to_structure_dictionary(self):
        """Add element specific keywords using the elementmap section
        """
        return

    def convert_int_to_atomidx(self, aint):
        """Convert an integer to an atom index.
            Args:
                aint <int>: integer
            Returns:
                atomidx <str>: atom index string
        """
        spad=16
        atomidx=hex(aint).zfill(spad)
        return atomidx

    def make_defects_dictionary(self, scaling_label=""):
        """Make a defect dictionary.
            Args: 
                scaling_label <str>: scaling label. Leave blank for no scaling.            
        """
        mydefdict=dict()
        if scaling == "":
            sdict = dict(self.startdict)
            mySE = self.startSE
        else:
            sdict = dict(self.scalingdicts[scaling_label])
            mySE = self.scalingSEs[scaling_label])
        #print input_options
        defect_dict=self.input_options.get_item('defects','defects')
        #print defect_dict
        dlabels=defect_dict.keys()
        for dlabel in dlabels:
            dsubkeys=defect_dict[dlabel].keys()
            for dsubkey in dsubkeys:
                if "subdefect_" in dsubkey:
                    dtype=defect_dict[dlabel][dsubkey]['type']
                    dcoords=defect_dict[dlabel][dsubkey]['coordinates']
                    if not (scaling_label == ""):
                        dcoords = mySE.get_scaled_coordinates(dcoords)
                    if dtype == "interstitial":
                        newdict=dict()
                        newdict['original_frac_coords']=dcoords
                        newdict['element']=defect_dict[dlabel][dsubkey]['symbol']
                        self.add_structure_dictionary_entry(sdict, newdict)
                    else:
                        didx=self.find_orig_frac_coord_in_structure_dictionary(sdict, dcoords)
                        if not dlabel in sdict.keys():
                            if dtype in ['substitution','antisite']:
                                sdict[didx]['element'] = defect_dict[dlabel][dsubkey]['symbol']
                            elif dtype == 'vacancy':
                                sdict.pop(didx,"None")
                        else:
                            raise MASTError(self.__class__.__name__,"Defect label %s already exists for atom index %s" % (dlabel, didx))
            self.write_structure_dictionary_file(sdict, "defect_%s%s_structure_index" % (otherlabels, dlabel))
        return sdict

    def find_orig_frac_coord_in_structure_dictionary(self, sdict, coord, tol=0.0001):
        """Find the atomic index of an original FRACTIONAL coordinate in the 
            structure dictionary.
            Args:
                sdict <dictionary>: structure dictioary
                coord <numpy array of float>: coordinate to find
                tol <float>: tolerance
            Returns:
                atomic index <hex string>: atomic index of match
                Returns None if no match is found
        """
        rtol=tol*100
        atomidxs=sdict.keys()
        matches=list()
        for atomidx in atomidxs:
            atom_ofc=sdict[atomidx]['original_frac_coords']
            if np.allclose(atom_ofc,coord,rtol,tol):
                matches.append(atomidx)
        if len(matches) > 1:
            raise MASTError(self.__class__.__name__,
                "Multiple matches found for coordinate %s" % coord)
        return matches[0]

        
    def add_structure_dictionary_entry(self, sdict, adict):
        """Add an entry to the structure dictionary, at the end.
            Args:
                sdict <dict>: Structure dictionary
                adict <dict>: Dictionary to add.
        """
        aidxs=sdict.keys()
        newidx=self.convert_int_to_atomidx(len(aidxs))
        if newidx in aidxs:
            raise MASTError(self.__class__.__name__, "Cannot overwrite atomic index %s" % newidx)
        sdict[newidx]=adict
        return True


    def write_structure_dictionary_file(self, sdict, sfilename):
        """Write a structure diciontary file from self.struc_dict
            Args:
                sdict <dict>: structure dictionary
                sfilename <str>: file name
        """
        satomindices=sdict.keys()
        fullheaders=list()
        for satomidx in satomindices:
            satomheaders=sdict[satomidx].keys()
            for satomheader in satomheaders:
                if not (satomheader in fullheaders):
                    fullheaders.append(satomheader)
        with open('%s' % sfilename, 'wb') as sfile:
            headerline='atomidx;'
            for fullheader in fullheaders:
                headerline = headerline + "%s;" % fullheader
            headerline = headerline + "\n"
            sfile.write(headerline)
            for satomidx in satomindices:
                atomline=""
                atomline="%s;" % satomidx
                satomkeys=sdict[satomidx].keys()
                for fullheader in fullheaders:
                    if fullheader in satomkeys:
                        atomline = atomline + "%s;" % sdict[satomidx][fullheader]
                    else:
                        atomline = atomline + ";"
                atomline = atomline + "\n"
                sfile.write(atomline)
        return
