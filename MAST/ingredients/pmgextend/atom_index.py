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
            'input_options': (InputOptions, None, 'Input options')
        }
        MASTObj.__init__(self, allowed_keys, **kwargs)            
        self.startSE = ""
        self.startdict = ""
        self.startdefects = ""
        self.startdefectphonons = ""
        self.startnebs = ""
        self.startnebphonons = ""
        self.scalingSEs = dict()
        self.scalingdicts = dict()
        self.scalingdefects = dict()
        self.scalingdefectphonons = dict()
        self.scalingnebs = dict()
        self.scalingnebphonons = dict()
        self.input_options = self.keywords['input_options']
        self.scaling = self.input_options.get_item('structure','scaling')
        if self.scaling == None:
            self.scaling = dict()
        self.atomcount=1
        return

    def set_up_initial_index(self):
        """Set up the initial index (first time).
        """
        self.get_structure_extensions()
        self.startdict = self.build_structure_dictionary(self.startSE)
        [self.startdefects, self.startdefectphonons] = self.make_defects_instruction_dictionary("")
        for scaling_label in self.scaling.keys():
            self.scalingdicts[scaling_label] = self.build_structure_dictionary(self.scalingSEs[scaling_label])
            [self.scalingdefects[scaling_label], self.scalingdefectphonons[scaling_label]] = self.make_defects_instruction_dictionary(scaling_label)
        print self.startdict
        print self.startdefects
        print self.startdefectphonons
        for scaling_label in self.scaling.keys():
            print self.scalingdicts[scaling_label]
            print self.scalingdefects[scaling_label]
            print self.scalingdefectphonons[scaling_label]
        return

    def get_structure_extensions(self):
        """Get starting structure extension object and scaled structure 
            extension objects
        """
        startstr=self.input_options.get_item('structure','structure')
        self.startSE = SE(struc_work1 = startstr)
        for scaling_label in self.scaling.keys():
            scaleSE = SE(struc_work1 = self.startSE.copy(), scaling_size=self.scaling[scaling_label][0])
            self.scalingSEs[scaling_label] = SE(struc_work1 = scaleSE.scale_structure())
        return 

    def build_structure_dictionary(self, mySE):
        """Build a structure dictionary.
        """
        sdict=dict()
        for site in mySE.keywords['struc_work1'].sites:
            skey=self.get_new_key() # starts at 1
            sdict[skey]=dict()
            sdict[skey]['original_frac_coords']=site.frac_coords
            sdict[skey]['element']=site.species_string
            sdict[skey]['specie']=site.specie
        return sdict

    def add_atom_specific_keywords_to_structure_dictionary(self):
        """Add atom specific keywords using the coordinates section
        """
        return
    
    def add_element_specific_keywords_to_structure_dictionary(self):
        """Add element specific keywords using the elementmap section
        """
        return

    def get_new_key(self):
        """Get a new key.
        """
        self.atomcount = self.atomcount + 1
        return self.convert_int_to_atomidx(self.atomcount)

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


    def make_defects_instruction_dictionary(self, scaling_label=""):
        """Make a defect instruction dictionary.
            Args: 
                scaling_label <str>: scaling label. Leave blank for no scaling.            
        """
        mydefdict=dict()
        myphondict=dict()
        if scaling_label == "":
            sdict = dict(self.startdict)
            mySE = self.startSE
        else:
            sdict = dict(self.scalingdicts[scaling_label])
            mySE = self.scalingSEs[scaling_label]
        defect_dict=self.input_options.get_item('defects','defects')
        dlabels=defect_dict.keys()
        for dlabel in dlabels:
            mydefdict[dlabel] = dict()
            mydefdict[dlabel]['add'] = dict()
            mydefdict[dlabel]['remove'] = dict()
            myphondict[dlabel] = dict()
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
                        
                        mydefdict[dlabel]['add'][self.get_new_key()]=newdict
                    else:
                        didx=self.find_orig_frac_coord_in_structure_dictionary(sdict, dcoords)
                        if dtype in ['substitution','antisite']:
                            mydefdict[dlabel]['remove'][didx] = 'remove'
                            newdict=dict()
                            newdict['original_frac_coords']=dcoords
                            newdict['element']=defect_dict[dlabel][dsubkey]['symbol']
                            mydefdict[dlabel]['add'][self.get_new_key()]=newdict
                        elif dtype == 'vacancy':
                            mydefdict[dlabel]['remove'][didx] = 'remove'
                if "phonon" in dsubkey:
                    for phonlabel in defect_dict[dlabel][dsubkey].keys():
                        pcoordsraw = defect_dict[dlabel][dsubkey][phonlabel]['phonon_center_site']
                        pthresh = defect_dict[dlabel][dsubkey][phonlabel]['threshold']
                        pcrad = defect_dict[dlabel][dsubkey][phonlabel]['phonon_center_radius']
                        pcoords = np.array(pcoordsraw.split(),'float')
                        
                        pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, pthresh+pcrad, True)
                        myphondict[dlabel][phonlabel] = list(pindices)
        return [mydefdict, myphondict]

    #'phonon': {'solute': {'phonon_center_site': '0.25 0.50 0.25', 'threshold': 0.1, 'phonon_center_radius': 0.5}

    def find_orig_frac_coord_in_structure_dictionary(self, sdict, coord, tol=0.0001, find_multiple=False):
        """Find the atomic index of an original FRACTIONAL coordinate in the 
            structure dictionary.
            Args:
                sdict <dictionary>: structure dictioary
                coord <numpy array of float>: coordinate to find
                tol <float>: tolerance
                find_multiple <boolean>: allow multiple matches. Default False.
            Returns:
                atomic index <hex string>: atomic index of match, 
                    if find_multiple is false
                list of atomic indices of matches, if find_multiple is true
                Returns None if no match is found
        """
        rtol=tol*100
        atomidxs=sdict.keys()
        matches=list()
        for atomidx in atomidxs:
            atom_ofc=sdict[atomidx]['original_frac_coords']
            if np.allclose(atom_ofc,coord,rtol,tol):
                matches.append(atomidx)
        if not find_multiple:
            if (len(matches) > 1):
                raise MASTError(self.__class__.__name__,
                    "Multiple matches found for coordinate %s: %s" % (coord, matches))
            if len(matches) == 1:
                return matches[0]
            else:
                return None
        else:
            if len(matches) >= 1:
                return matches
            else:
                return None

        
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
