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
        self.allatoms = ""
        return

    def set_up_initial_index(self):
        """Set up the initial index (first time).
        """
        self.get_structure_extensions()
        self.startdict = self.build_structure_dictionary(self.startSE)
        [self.startdefects, self.startdefectphonons] = self.make_defects_instruction_dictionary("")
        [self.startnebs, self.startnebphonons] = self.make_neb_instruction_dictionary("")
        for scaling_label in self.scaling.keys():
            self.scalingdicts[scaling_label] = self.build_structure_dictionary(self.scalingSEs[scaling_label])
            [self.scalingdefects[scaling_label], self.scalingdefectphonons[scaling_label]] = self.make_defects_instruction_dictionary(scaling_label)
            [self.scalingnebs[scaling_label], self.scalingnebphonons[scaling_label]] = self.make_neb_instruction_dictionary(scaling_label)
        print self.startdict
        print self.startdefects
        print self.startdefectphonons
        print self.startnebs
        print self.startnebphonons
        for scaling_label in self.scaling.keys():
            print self.scalingdicts[scaling_label]
            print self.scalingdefects[scaling_label]
            print self.scalingdefectphonons[scaling_label]
            print self.scalingnebs[scaling_label]
            print self.scalingnebphonons[scaling_label]
        self.allatoms = self.combine_structure_dictionaries()
        print self.allatoms
        print self.input_options
        self.make_structure_index_directory()
        self.make_manifest_files()
        return

    def make_structure_index_directory(self):
        """Make structure index directory
        """
        if os.path.isdir("structure_index_files"):
            import time
            time.sleep(1)
        if os.path.isdir("structure_index_files"):
            raise MASTError(self.__class__.__name__, "Structure index directory already exists!")
        os.mkdir("structure_index_files")
        self.write_atom_index_files()
        return

    def write_atom_index_files(self):
        """Write atomic index files
        """
        for akey in self.allatoms:
            aname = os.path.join("structure_index_files","atom_index_%s" % akey)
            afile = open(aname,'wb')
            afile.write("%20s:%20s\n" % ("name",akey))
            for key, value in self.allatoms[akey].iteritems():
                afile.write("%20s:%20s\n" % (key, value))
            afile.close()
        return

    def read_atom_index_file(self, filename):
        """Read an atom index file.
            Args:
                filename <str>: File name to read
        """
        adict=dict()
        afile = open(filename, 'rb')
        alines = afile.readlines()
        afile.close()
        for aline in alines:
            asplit = aline.split(":",1)
            akey = asplit[0].strip()
            aval = asplit[1].strip()
            adict[akey] = aval
        return adict

    def make_manifest_file(self, scaling_label="", defect_label="", neb_label="", phonon_label=""):
        """Make manifest file.
            Args:
                scaling_label <str>: Scaling label
                defect_label <str>: Defect label. NEBs should ALSO contain
                            a defect_label to specify which endpoint manifest
                            is being created
                neb_label <str>: NEB label
                phonon_label <str> phonon label
        """
        fname="manifest:%s:%s:%s:%s" % (scaling_label, defect_label, neb_label, phonon_label)
        fname = os.path.join("structure_index_files", fname)
        if scaling_label == "":
            sdict = dict(self.startdict)
            defdict = dict(self.startdefects)
            defphondict = dict(self.startdefectphonons)
            nebdict = dict(self.startnebs)
            nebphondict = dict(self.startnebphonons)
        else:
            sdict = dict(self.scalingdicts[scaling_label])
            defdict = dict(self.scalingdefects[scaling_label])
            defphondict = dict(self.scalingdefectphonons[scaling_label])
            nebdict = dict(self.scalingnebs[scaling_label])
            nebphondict = dict(self.scalingnebphonons[scaling_label])
        alist = list(sdict.keys()) #typically these are already in element order
        if not (defect_label == ""):
            for addme in defdict[defect_label]['add'].keys():
                alist.append(addme)
            for removeme in defdict[defect_label]['remove'].keys():
                alist.remove(removeme)
        if not (neb_label == ""):
            matchlist = list(nebdict[neb_label]['match'])
            nebsplit = neb_label.split("-")
            if defect_label in nebsplit[0]:
                whichep = 0
            else:
                whichep = 1
            for matchline in matchlist:
                alist.remove(matchline[whichep]) #remove from middle
                alist.append(matchline[whichep]) #add to the bottom
        self.write_manifest_file(alist, fname)
        return

    def make_manifest_files(self):
        """Make manifest files.
        """
        self.make_manifest_file()
        for defect_label in self.startdefects.keys():
            self.make_manifest_file("",defect_label)
            for scaling_label in self.scaling.keys():
                self.make_manifest_file(scaling_label,defect_label)
            for neb_label in self.startnebs.keys():
                self.make_manifest_file("",defect_label,neb_label)
                for scaling_label in self.scaling.keys():
                    self.make_manifest_file(scaling_label,defect_label,neb_label)
        return

    def get_atoms_for_phonons(self):
        """
        """
        return

    def OLD_set_structure_from_inputs(self, input_options):
        """Make a pymatgen structure and update the
            structure key.
            Args:
                input_options <InputOptions>
        """
        strposfile = input_options.get_item('structure','posfile')
        if strposfile is None:
            iopscoords=input_options.get_item('structure','coordinates')
            iopslatt=input_options.get_item('structure','lattice')
            iopsatoms=input_options.get_item('structure','atom_list')
            iopsctype=input_options.get_item('structure','coord_type')
            structure = MAST2Structure(lattice=iopslatt,
                coordinates=iopscoords, atom_list=iopsatoms,
                coord_type=iopsctype)
        elif ('poscar' in strposfile.lower()):
            from pymatgen.io.vaspio import Poscar
            structure = Poscar.from_file(strposfile).structure
        elif ('cif' in strposfile.lower()):
            from pymatgen.io.cifio import CifParser
            structure = CifParser(strposfile).get_structures()[0]
        else:
            error = 'Cannot build structure from file %s' % strposfile
            raise MASTError(self.__class__.__name__, error)
        input_options.update_item('structure','structure',structure)
        if not input_options.get_item('structure','use_structure_index'):
            pass
        else:
            self.do_structure_indexing(input_options)
        return

    def write_manifest_file(self, aidxlist, fname):
        """Make a manifest file.
            Args:
                aidxlist <list of str>: List of atom indices
                fname <str>: File name
        """
        myfile = MASTFile()
        myfile.data = aidxlist
        myfile.to_file(fname)
        return

    
    def combine_structure_dictionaries(self):
        """Combine structure dictionaries into single comprehensive dictionary.
        """
        largedict=dict()
        for skey in self.startdict.keys():
            largedict[skey]=dict(self.startdict[skey])
        for dkey in self.startdefects.keys():
            for addkey in self.startdefects[dkey]['add'].keys():
                largedict[addkey]=dict(self.startdefects[dkey]['add'][addkey])
        for scaling_label in self.scaling.keys():
            for skey in self.scalingdicts[scaling_label].keys():
                largedict[skey]= dict(self.scalingdicts[scaling_label][skey])
            for dkey in self.scalingdefects[scaling_label].keys():
                for addkey in self.scalingdefects[scaling_label][dkey]['add'].keys():
                    largedict[addkey]=dict(self.scalingdefects[scaling_label][dkey]['add'][addkey])
        return largedict

    def get_structure_extensions(self):
        """Get starting structure extension object and scaled structure 
            extension objects
        """
        startstr=self.input_options.get_item('structure','structure')
        self.startSE = SE(struc_work1 = startstr)
        for scaling_label in self.scaling.keys():
            scaleSE = SE(struc_work1 = self.startSE.keywords['struc_work1'].copy(), scaling_size=self.scaling[scaling_label][0])
            self.scalingSEs[scaling_label] = scaleSE
        return 

    def build_structure_dictionary(self, mySE):
        """Build a structure dictionary.
        """
        sdict=dict()
        if 'scaling_label' in mySE.keywords.keys():
            mystruc = mySE.scale_structure()
        else:
            mystruc = mySE.keywords['struc_work1']
        for site in mystruc:
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
        """Make a defect instruction dictionary and a defects phonon dictionary.
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
        if defect_dict == None:
            return [dict(), dict()]
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
    
    #'mast_neb_settings': {'images': 1, 'phonon': {'movingsolvent': {'phonon_center_site': '0.375 0.5 0.375', 'threshold': 0.1, 'phonon_center_radius': 0.5}}, 'lines': [['Al', ' 0.25 0.5 0.25', ' 0.5 0.5 0.5']]}, 'mast_ppn': '1', 'ismear': '1', 'nebs': {'1nn-solute': {'images': 1, 'phonon': {'movingsolute': {'phonon_center_site': '0.375 0.500 0.375', 'threshold': 0.1, 'phonon_center_radius': 0.5}}, 'lines': [['Mg', ' 0.25 0.50 0.25', ' 0.5 0.5 0.5']]}, 'pureinit-purefin': {'images': 1, 'phonon': {'movingsolvent': {'phonon_center_site': '0.375 0.5 0.375', 'threshold': 0.1, 'phonon_center_radius': 0.5}}, 'lines': [['Al', ' 0.25 0.5 0.25', ' 0.5 0.5 0.5']]}}
    def make_neb_instruction_dictionary(self, scaling_label=""):
        """Make a neb instruction dictionary and a neb phonon dictionary.
            Args: 
                scaling_label <str>: scaling label. Leave blank for no scaling.            
        """
        mynebdict=dict()
        myphondict=dict()
        if scaling_label == "":
            sdict = dict(self.startdict)
            mySE = self.startSE
        else:
            sdict = dict(self.scalingdicts[scaling_label])
            mySE = self.scalingSEs[scaling_label]
        neb_dict=self.input_options.get_item('neb','nebs')
        if neb_dict == None:
            return [dict(), dict()]
        nlabels=neb_dict.keys()
        for nlabel in nlabels:
            mynebdict[nlabel] = dict()
            mynebdict[nlabel]['match'] = list()
            myphondict[nlabel] = dict()
            nsubkeys=neb_dict[nlabel].keys()
            for nsubkey in nsubkeys:
                if "lines" in nsubkey:
                    nlines = list(neb_dict[nlabel][nsubkey])
                    for nline in nlines:
                        ncoord1 = np.array(nline[1].split(), 'float')
                        ncoord2 = np.array(nline[2].split(), 'float')
                        nidx1 = self.find_orig_frac_coord_in_structure_dictionary(sdict, ncoord1)
                        nidx2 = self.find_orig_frac_coord_in_structure_dictionary(sdict, ncoord2)
                        mynebdict[nlabel]['match'].append([nidx1, nidx2])
                if "phonon" in nsubkey:
                    for phonlabel in neb_dict[nlabel][nsubkey].keys():
                        pcoordsraw = neb_dict[nlabel][nsubkey][phonlabel]['phonon_center_site']
                        pthresh = neb_dict[nlabel][nsubkey][phonlabel]['threshold']
                        pcrad = neb_dict[nlabel][nsubkey][phonlabel]['phonon_center_radius']
                        pcoords = np.array(pcoordsraw.split(),'float')
                        
                        pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, pthresh+pcrad, True)
                        myphondict[nlabel][phonlabel] = list(pindices)
        return [mynebdict, myphondict]

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
