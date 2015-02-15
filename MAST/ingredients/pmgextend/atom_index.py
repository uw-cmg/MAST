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
from MAST.utility import Metadata
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
        self.sdir = "structure_index_files"
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
        self.startstr = self.input_options.get_item('structure','structure')
        self.atomcount=1
        self.allatoms = ""
        return
    
    def make_structure_index_directory(self):
        """Make structure index directory
        """
        if os.path.isdir(self.sdir):
            import time
            time.sleep(1)
        if os.path.isdir(self.sdir):
            raise MASTError(self.__class__.__name__, "Structure index directory already exists!")
        os.mkdir(self.sdir)
        return

    def write_undefected_atom_indices(self):
        """Write undefected atom indices, including scaled indices.
            Also write an undefected manifest file.
        """
        scales = self.scaling.keys()
        scales.append("")
        for scaling_label in scales:
            if scaling_label == "":
                mySE=SE(struc_work1=self.startstr.copy())
                mystruc=mySE.keywords['struc_work1']
            else:
                mySE=SE(struc_work1=self.startstr.copy(), scaling_size=self.scaling[scaling_label][0])
                mystruc=mySE.scale_structure()
            alist=list()
            manname=os.path.join(self.sdir,"manifest_undefected_%s" % scaling_label)
            for site in mystruc:
                akey=self.get_new_key()
                aname="atomindex_%s" % akey
                aname = os.path.join(self.sdir, aname)
                ameta = Metadata(metafile=aname)
                ameta.write_data("atomindex",akey)
                ameta.write_data("original_frac_coords", site.frac_coords)
                ameta.write_data("element", site.species_string)
                ameta.write_data("scaling_label", scaling_label)
                alist.append(akey)
            self.write_manifest_file(alist,manname)
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
    
    def write_manifest_file(self, aidxlist, fname):
        """Make a manifest file.
            Args:
                aidxlist <list of str>: List of atom indices
                fname <str>: File name
        """
        myfile = open(fname, 'wb')
        for aidx in aidxlist:
            myfile.write("%s\n" % aidx)
        myfile.close()
        return

    def write_defected_atom_indices(self):
        """Write any additional defect atom indices and make manifests.
        """
        defect_dict=self.input_options.get_item('defects','defects')
        if defect_dict == None:
            return None
        dlabels=defect_dict.keys()
        
        scales = self.scaling.keys()
        scales.append("")
        for scaling_label in scales:
            alist=list(self.read_manifest_file("%s/manifest_undefected_%s" % (self.sdir, scaling_label)))
            if scaling_label == "":
                mySE=SE(struc_work1=self.startstr.copy())
            else:
                mySE=SE(struc_work1=self.startstr.copy(), scaling_size=self.scaling[scaling_label][0])
            for dlabel in dlabels:
                dlist = list(alist)
                manname=os.path.join(self.sdir,"manifest_defect_%s_%s" % (dlabel, scaling_label))
                dsubkeys=defect_dict[dlabel].keys()
                for dsubkey in dsubkeys:
                    if "subdefect_" in dsubkey:
                        dtype=defect_dict[dlabel][dsubkey]['type']
                        dcoords=defect_dict[dlabel][dsubkey]['coordinates']
                        delement=defect_dict[dlabel][dsubkey]['symbol']
                        if not (scaling_label == ""):
                            dcoords = mySE.get_scaled_coordinates(dcoords)
                        if dtype == "interstitial":
                            didx=self.find_orig_frac_coord_in_atom_indices(dcoords, delement, scaling_label, False, 0.001)
                            if didx == None:
                                akey=self.get_new_key()
                                aname="atomindex_%s" % akey
                                aname = os.path.join(self.sdir, aname)
                                ameta = Metadata(metafile=aname)
                                ameta.write_data("atomindex",akey)
                                ameta.write_data("original_frac_coords", dcoords)
                                ameta.write_data("element", delement)
                                ameta.write_data("scaling_label", scaling_label)
                                dlist.append(akey)
                            else:
                                dlist.append(didx)
                        elif dtype == "vacancy":
                            didx=self.find_orig_frac_coord_in_atom_indices(dcoords, delement, scaling_label, False, 0.001)
                            dlist.remove(didx)
                        elif dtype in ["substitution","antisite"]:
                            didxlist=self.find_orig_frac_coord_in_atom_indices(dcoords, "", scaling_label, True, 0.001) #leave element empty; just search coords
                            for didx in didxlist:
                                dmeta = Metadata(metafile="%s/atom_index_%s" % (self.sdir, didx))
                                dmetaelem = dmeta.read_data("element")
                                if not (delement == dmetaelem):
                                    if didx in dlist:
                                        dlist.remove(didx)
                            didxsub=self.find_orig_frac_coord_in_atom_indices(dcoords, delement, scaling_label, False, 0.001) #leave element empty; just search coords
                            if didxsub == None:
                                akey=self.get_new_key()
                                aname="atomindex_%s" % akey
                                aname = os.path.join(self.sdir, aname)
                                ameta = Metadata(metafile=aname)
                                ameta.write_data("atomindex",akey)
                                ameta.write_data("original_frac_coords", dcoords)
                                ameta.write_data("element", delement) #sub element here
                                ameta.write_data("scaling_label", scaling_label)
                                dlist.append(akey)
                            else:
                                dlist.append(didxsub)
                self.write_manifest_file(dlist, manname)
        return 
    
    def read_manifest_file(self, filename):
        """Read a manifest file.
        """
        mlist=list()
        mfile = open(filename, 'rb')
        mlines = mfiles.readlines()
        mfile.close()
        for mline in mlines:
            mline = mline.strip()
            mlist.append(mline)
        return mlist
    
    def find_orig_frac_coord_in_atom_indices(self, coord, element="", scaling_label="", find_multiple=False, tol=0.0001):
        """Find the atomic index of an original FRACTIONAL coordinate in the 
            structure dictionary.
            Args:
                coord <numpy array of float>: coordinate to find
                element <str>: element symbol to match
                scaling_label <str>: scaling label ("" for no scaling)
                find_multiple <boolean>: allow multiple matches. Default False.
                tol <float>: tolerance
            Returns:
                atomic index <hex string>: atomic index of match, 
                    if find_multiple is false
                list of atomic indices of matches, if find_multiple is true
                Returns None if no match is found
        """
        import glob
        matchstring = "%s/atom_index_*" % self.sdir
        idxnames = glob.glob(matchstring)
        rtol=tol*100
        coord_matches=list()
        elem_matches=list()
        scaling_matches=list()
        for aname in idxnames:
            ameta=Metadata(metafile=aname)
            aidx=ameta.read_data("atomindex")
            atom_ofc=ameta.read_data("original_frac_coords")
            atom_ofc_arr=np.array(atom_ofc[1:-1].split(),'float')
            if np.allclose(atom_ofc_arr,coord,rtol,tol):
                coord_matches.append(aidx)
        if element == "":
            elem_matches = list(coord_matches)
        else:
            for aidx in coord_matches:
                ameta=Metadata(metafile="%s/atom_index_%s" % (self.sdir, aidx))
                atom_elem=ameta.read_data("element")
                if (element == atom_elem):
                    elem_matches.append(aidx)
        if scaling_label == "":
            scaling_matches = list(elem_matches)
        else:
            for aidx in elem_matches:
                ameta=Metadata(metafile="%s/atom_index_%s" % (self.sdir, aidx))
                ascale=ameta.read_data("scaling_label")
                if (scaling_label == ascale):
                    scaling_matches.append(aidx)
        allmatches = list(scaling_matches)
        if len(allmatches) == 0:
            return None
        if len(allmatches) > 1:
            if not find_multiple:
                raise MASTError(self.__class__.__name__,
                    "Multiple matches found for coordinate %s: %s" % (coord, allmatches))
            else:
                return allmatches
        if len(allmatches) == 1:
            if not find_multiple:
                return allmatches[0]
            else:
                return allmatches
       return None

    def write_defected_phonon_sd_manifests(self):
        """Write defected phonon structure dynamics manifests.
        """
        defect_dict=self.input_options.get_item('defects','defects')
        if defect_dict == None:
            return None
        dlabels=defect_dict.keys()
        
        scales = self.scaling.keys()
        scales.append("")
        for scaling_label in scales:
            if scaling_label == "":
                mySE=SE(struc_work1=self.startstr.copy())
            else:
                mySE=SE(struc_work1=self.startstr.copy(), scaling_size=self.scaling[scaling_label][0])
            for dlabel in dlabels:
                pdict=dict(defect_dict[dlabel]["phonon"])
                for phonon_label in pdict.keys():
                    pcoordsraw = pdict[phonon_label]['phonon_center_site']
                    pthresh = pdict[phonon_label]['threshold']
                    pcrad = pdict[phonon_label]['phonon_center_radius']
                    pcoords = np.array(pcoordsraw.split(),'float')
                    if not (scaling_label == ""):
                        pcoords = mySE.get_scaled_coordinates(pcoords)
                     
                    #pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, pthresh+pcrad, True)
                    pindices = self.find_orig_frac_coord_in_atom_indices(pcoords,"",scaling_label,True,0.001+pcrad)
                    manname=os.path.join(self.sdir,"manifest_phonon_sd_%s_%s" % (dlabel, phonon_label, scaling_label))
                    self.write_manifest_file(pindices, manname) 
        return 

    def write_neb_endpoint_manifests():
        """Make NEB endpoint manifests.
        """
        neb_dict=self.input_options.get_item('neb','nebs')
        if neb_dict == None:
            return None
        nlabels=neb_dict.keys()
        
        scales = self.scaling.keys()
        scales.append("")
        for scaling_label in scales:
            if scaling_label == "":
                mySE=SE(struc_work1=self.startstr.copy())
            else:
                mySE=SE(struc_work1=self.startstr.copy(), scaling_size=self.scaling[scaling_label][0])
            for nlabel in nlabels:
                def1 = nlabel.split("-")[0].strip()
                def2 = nlabel.split("-")[1].strip()
                manname1=os.path.join(self.sdir,"manifest_neb_%s_%s_%s" % (nlabel, def1, scaling_label))
                manname2=os.path.join(self.sdir,"manifest_neb_%s_%s_%s" % (nlabel, def2, scaling_label))
                mlist1=list(self.read_manifest_file("%s/manifest_defect_%s_%s" % (self.sdir, def1, scaling_label)))
                mlist2=list(self.read_manifest_file("%s/manifest_defect_%s_%s" % (self.sdir, def2, scaling_label)))
                maddtoend1=list()
                maddtoend2=list()
                nlines=list(neb_dict[nlabel]["lines"])
                for nline in nlines:
                    ncoord1 = np.array(nline[1].split(), 'float')
                    ncoord2 = np.array(nline[2].split(), 'float')
                    if not (scaling_label == ""):
                        ncoord1 = mySE.get_scaled_coordinates(ncoord1)
                        ncoord2 = mySE.get_scaled_coordinates(ncoord2)
                    nelem = nline[0]
                    nidx1 = self.find_orig_frac_coord_in_atom_indices(ncoord1, nelem, scaling_label, False, 0.001)
                    nidx2 = self.find_orig_frac_coord_in_atom_indices(ncoord1, nelem, scaling_label, False, 0.001)
                    mlist1.remove(nidx1)
                    maddtoend1.append(nidx1) #resort matches to the bottom
                    mlist2.remove(nidx2)
                    maddtoend2.append(nidx2)
                for midx in [0:len(mlist1)]:
                    if not (mlist1[midx] == mlist2[midx]):
                        raise MASTError("NEB %s truncated manifests do not match: %s, %s" % (nlabel, mlist1, mlist2))
                mlist1.extend(maddtoend1)
                mlist2.extend(maddtoend2)
                self.write_manifest_file(mlist1, manname1)
                self.write_manifest_file(mlist2, manname2)
        return
    def write_neb_phonon_sd_manifests():
        """Make NEB phonon manifests.
        """
        neb_dict=self.input_options.get_item('neb','nebs')
        if neb_dict == None:
            return None
        nlabels=neb_dict.keys()
        
        scales = self.scaling.keys()
        scales.append("")
        for scaling_label in scales:
            if scaling_label == "":
                mySE=SE(struc_work1=self.startstr.copy())
            else:
                mySE=SE(struc_work1=self.startstr.copy(), scaling_size=self.scaling[scaling_label][0])
            for nlabel in nlabels:
                def1 = nlabel.split("-")[0].strip()
                def2 = nlabel.split("-")[1].strip()
                manname1=os.path.join(self.sdir,"manifest_neb_%s_%s_%s" % (nlabel, def1, scaling_label))
                manname2=os.path.join(self.sdir,"manifest_neb_%s_%s_%s" % (nlabel, def2, scaling_label))
                mlist1=list(self.read_manifest_file("%s/manifest_defect_%s_%s" % (self.sdir, def1, scaling_label)))
                mlist2=list(self.read_manifest_file("%s/manifest_defect_%s_%s" % (self.sdir, def2, scaling_label)))
                nlines=list(neb_dict[nlabel]["lines"])
                for nline in nlines:
                    ncoord1 = np.array(nline[1].split(), 'float')
                    ncoord2 = np.array(nline[2].split(), 'float')
                    if not (scaling_label == ""):
                        ncoord1 = mySE.get_scaled_coordinates(ncoord1)
                        ncoord2 = mySE.get_scaled_coordinates(ncoord2)
                    nelem = nline[0]
                    nidx1 = self.find_orig_frac_coord_in_atom_indices(ncoord1, nelem, scaling_label, False, 0.001)
                    nidx2 = self.find_orig_frac_coord_in_atom_indices(ncoord1, nelem, scaling_label, False, 0.001)
                    mlist1.remove(nidx1)
                    mlist1.append(nidx1) #resort to the bottom
                    mlist2.remove(nidx2)
                    mlist2.append(nidx2)
                self.write_manifest_file(mlist1, manname1)
                self.write_manifest_file(mlist2, manname2)
        return



    def set_up_initial_index(self):
        """Set up the initial index (first time).
        """
        self.make_structure_index_directory()
        self.write_undefected_atom_indices()
        self.write_defected_atom_indices()
        self.write_defected_phonon_sd_manifests()
        self.write_neb_endpoint_manifests()
        
        ##
        self.get_structure_extensions()
        self.startdict = self.build_structure_dictionary(self.startSE)
        self.startdefects = self.make_defects_instruction_dictionary("")
        for scaling_label in self.scaling.keys():
            self.scalingdicts[scaling_label] = self.build_structure_dictionary(self.scalingSEs[scaling_label])
            self.scalingdefects[scaling_label] = self.make_defects_instruction_dictionary(scaling_label)
        self.allatoms = self.combine_structure_dictionaries()
        self.startdefectphonons = self.make_defects_phonon_dictionary("")
        [self.startnebs, self.startnebphonons] = self.make_neb_instruction_dictionary("")
        for scaling_label in self.scaling.keys():
            self.scalingdefectphonons[scaling_label] = self.make_defects_phonon_dictionary(scaling_label)
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
        print self.allatoms
        print self.input_options
        return


    def write_initial_manifests(self):
        """Write initial structure manifests.
        """
        self.write_single_initial_manifest()
        for scaling_label in self.scaling.keys():
            self.write_single_initial_manifest(scaling_label)
        return

    def write_single_initial_manifest(self, scaling_label=""):
        """Write a single manifest file
        """
        fname="initial_manifest_%s" % scaling_label
        fname = os.path.join("structure_index_files", fname)
        if scaling_label == "":
            sdict = dict(self.startdict)
        else:
            sdict = dict(self.scalingdicts[scaling_label])
        alist = list(sdict.keys()) #typically these are already in element order
        alist.sort() #sorts IN PLACE
        self.write_manifest_file(alist, fname)
        return

    def write_defect_instructions(self):
        """
        """
        for defect_label in self.startdefects.keys():
            self.write_single_defect_instruction(defect_label)
            for scaling_label in self.scaling.keys():
                self.write_single_defect_instruction(defect_label, scaling_label)
        return

    def write_single_defect_instruction(self, defect_label, scaling_label=""):
        """
            Args:
                defect_label <str>
                scaling_label <str>
        """
        dname="defect_instructions_%s_%s" % (defect_label, scaling_label)
        dname = os.path.join("structure_index_files", dname)
        if scaling_label == "":
            ddict = dict(self.startdefects[defect_label])
            pdict = dict(self.startdefectphonons[defect_label])
        else:
            ddict = dict(self.scalingdefects[scaling_label][defect_label])
            pdict = dict(self.scalingdefectphonons[scaling_label][defect_label])
        dfile = open(dname, 'wb')
        for arkey in ddict.keys():
            if arkey == 'replace':
                for aidxkey in ddict[arkey].keys():
                    for repkey in ddict[arkey][aidxkey].keys():
                        dfile.write("%s:%s:with:%s\n" % (arkey, aidxkey, repkey))
            else:
                for aidxkey in ddict[arkey].keys():
                    dfile.write("%s:%s\n" % (arkey, aidxkey))
        dfile.close()
        for pkey in pdict.keys():
            dpname="phonon_instructions_%s_%s_%s" % (defect_label, pkey, scaling_label)
            dpname = os.path.join("structure_index_files", dpname)
            dpfile = open(dpname, 'wb')
            for pentry in pdict[pkey]: #this is a list
                dpfile.write("%s\n" % (pentry))
            dpfile.close()
        return

    def write_neb_instructions(self):
        """
        """
        for neb_label in self.startnebs.keys():
            self.write_single_neb_instruction(neb_label)
            for scaling_label in self.scaling.keys():
                self.write_single_neb_instruction(neb_label, scaling_label)
        return

    def write_single_neb_instruction(self, neb_label, scaling_label=""):
        """
            Args:
                neb_label <str>
                scaling_label <str>
        """
        nname="neb_instructions_%s_%s" % (neb_label, scaling_label)
        nname = os.path.join("structure_index_files", nname)
        if scaling_label == "":
            ndict = dict(self.startnebs[neb_label])
            npdict = dict(self.startnebphonons[neb_label])
        else:
            ndict = dict(self.scalingnebs[scaling_label][neb_label])
            npdict = dict(self.scalingnebphonons[scaling_label][neb_label])
        nfile = open(nname, 'wb')
        for matchkey in ndict.keys():
            for matchentry in ndict[matchkey]: # this is a list
                nfile.write("%s:%s\n" % (matchentry[0],matchentry[1]))
        nfile.close()
        for pkey in npdict.keys():
            npname="phonon_instructions_%s_%s_%s" % (neb_label, pkey, scaling_label)
            npname = os.path.join("structure_index_files", npname)
            npfile = open(npname, 'wb')
            for pentry in npdict[pkey]: #this is a list
                npfile.write("%s\n" % (pentry))
            npfile.close()
        return
    


    def write_atom_index_files(self):
        """Write atomic index files
        """
        for akey in self.allatoms.keys():
            self.write_atom_index_file(akey)
        return

    def write_atom_index_file(self, aidx):
        """Write a single atom index file.
            Args:
                aidx <str>: atom index
        """
        aname = os.path.join("structure_index_files","atom_index_%s" % aidx)
        ameta = Metadata(metafile=aname)
        ameta.write_data("atomindex", aidx)
        for key, value in self.allatoms[aidx].iteritems():
            ameta.write_data(key, value)
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
            asplit = aline.split("=",1)
            akey = asplit[0].strip()
            aval = asplit[1].strip()
            adict[akey] = aval
        return adict

    
    def error_check_print_list(self, mylist, mylabel):
        """
        """
        print "LIST: ", mylabel
        for myentry in mylist:
            print myentry
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


    
    def combine_structure_dictionaries(self):
        """Combine structure dictionaries into single comprehensive dictionary.
        """
        largedict=dict()
        for skey in self.startdict.keys():
            largedict[skey]=dict(self.startdict[skey])
        for dkey in self.startdefects.keys():
            for addkey in self.startdefects[dkey]['add'].keys():
                largedict[addkey]=dict(self.startdefects[dkey]['add'][addkey])
            for repkey in self.startdefects[dkey]['replace'].keys():
                for repwithkey in self.startdefects[dkey]['replace'][repkey].keys():
                    largedict[repwithkey]=dict(self.startdefects[dkey]['replace'][repkey][repwithkey])
        for scaling_label in self.scaling.keys():
            for skey in self.scalingdicts[scaling_label].keys():
                largedict[skey]= dict(self.scalingdicts[scaling_label][skey])
            for dkey in self.scalingdefects[scaling_label].keys():
                for addkey in self.scalingdefects[scaling_label][dkey]['add'].keys():
                    largedict[addkey]=dict(self.scalingdefects[scaling_label][dkey]['add'][addkey])
                for repkey in self.scalingdefects[scaling_label][dkey]['replace'].keys():
                    for repwithkey in self.scalingdefects[scaling_label][dkey]['replace'][repkey].keys():
                        largedict[repwithkey]=dict(self.scalingdefects[scaling_label][dkey]['replace'][repkey][repwithkey])
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
        if mySE.keywords['scaling_size'] == None:
            mystruc = mySE.keywords['struc_work1']
        else:
            mystruc = mySE.scale_structure()
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



    def make_defects_instruction_dictionary(self, scaling_label=""):
        """Make a defect instruction dictionary 
            Args: 
                scaling_label <str>: scaling label. Leave blank for no scaling.            
        """
        mydefdict=dict()
        if scaling_label == "":
            sdict = dict(self.startdict)
            mySE = self.startSE
        else:
            sdict = dict(self.scalingdicts[scaling_label])
            mySE = self.scalingSEs[scaling_label]
        defect_dict=self.input_options.get_item('defects','defects')
        if defect_dict == None:
            return None
        dlabels=defect_dict.keys()
        for dlabel in dlabels:
            mydefdict[dlabel] = dict()
            mydefdict[dlabel]['add'] = dict()
            mydefdict[dlabel]['replace'] = dict()
            mydefdict[dlabel]['remove'] = dict()
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
                            newdict=dict()
                            newdict['original_frac_coords']=dcoords
                            newdict['element']=defect_dict[dlabel][dsubkey]['symbol']
                            mydefdict[dlabel]['replace'][didx]=dict()
                            mydefdict[dlabel]['replace'][didx][self.get_new_key()]=newdict
                        elif dtype == 'vacancy':
                            mydefdict[dlabel]['remove'][didx] = 'remove'
        return mydefdict
    
    def make_defects_phonon_dictionary(self, scaling_label=""):
        """
        """
        myphondict=dict()
        sdict = dict(self.allatoms)
        if scaling_label == "":
            mySE = self.startSE
        else:
            mySE = self.scalingSEs[scaling_label]
        defect_dict=self.input_options.get_item('defects','defects')
        if defect_dict == None:
            return dict()
        dlabels=defect_dict.keys()
        for dlabel in dlabels:
            myphondict[dlabel] = dict()
            dsubkeys=defect_dict[dlabel].keys()
            for dsubkey in dsubkeys:
                if "phonon" in dsubkey:
                    for phonlabel in defect_dict[dlabel][dsubkey].keys():
                        pcoordsraw = defect_dict[dlabel][dsubkey][phonlabel]['phonon_center_site']
                        pthresh = defect_dict[dlabel][dsubkey][phonlabel]['threshold']
                        pcrad = defect_dict[dlabel][dsubkey][phonlabel]['phonon_center_radius']
                        pcoords = np.array(pcoordsraw.split(),'float')
                        if not (scaling_label == ""):
                            pcoords = mySE.get_scaled_coordinates(pcoords)
                         
                        #pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, pthresh+pcrad, True)
                        pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, 0.001+pcrad, True)
                        myphondict[dlabel][phonlabel] = list(pindices)
        return myphondict

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
            sdict = dict(self.allatoms)
            mySE = self.startSE
        else:
            sdict = dict(self.allatoms)
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
                        if not (scaling_label == ""):
                            ncoord1 = mySE.get_scaled_coordinates(ncoord1)
                            ncoord2 = mySE.get_scaled_coordinates(ncoord2)
                        nelem = nline[0]
                        nidx1 = self.find_orig_frac_coord_in_structure_dictionary(sdict, ncoord1, 0.001, False, nelem) 
                        nidx2 = self.find_orig_frac_coord_in_structure_dictionary(sdict, ncoord2, 0.001, False, nelem)
                        mynebdict[nlabel]['match'].append([nidx1, nidx2])
                if "phonon" in nsubkey:
                    for phonlabel in neb_dict[nlabel][nsubkey].keys():
                        pcoordsraw = neb_dict[nlabel][nsubkey][phonlabel]['phonon_center_site']
                        pthresh = neb_dict[nlabel][nsubkey][phonlabel]['threshold']
                        pcrad = neb_dict[nlabel][nsubkey][phonlabel]['phonon_center_radius']
                        pcoords = np.array(pcoordsraw.split(),'float')
                        if not (scaling_label == ""):
                            pcoords = mySE.get_scaled_coordinates(pcoords)
                        
                        #pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, pthresh+pcrad, True)
                        pindices = self.find_orig_frac_coord_in_structure_dictionary(sdict, pcoords, 0.001+pcrad, True)
                        myphondict[nlabel][phonlabel] = list(pindices)
        return [mynebdict, myphondict]

    def find_orig_frac_coord_in_structure_dictionary(self, sdict, coord, tol=0.0001, find_multiple=False, element=""):
        """Find the atomic index of an original FRACTIONAL coordinate in the 
            structure dictionary.
            Args:
                sdict <dictionary>: structure dictioary
                coord <numpy array of float>: coordinate to find
                tol <float>: tolerance
                find_multiple <boolean>: allow multiple matches. Default False.
                element <str>: element symbol to match
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
                if element == "":
                    matches.append(atomidx)
                else:
                    if element == sdict[atomidx]['element']:
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

        


