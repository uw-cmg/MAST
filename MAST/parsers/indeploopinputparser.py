##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility.mastfile import MASTFile

ALLOWED_KEYS = {\
                 'inputfile'    : (str, 'mast.inp', 'Input file name'),\
               }

class IndepLoopInputParser(MASTObj):

    """Scans an input file for "indeploop" keyword and copies it into
        many input files.
        Attributes:
            self.indeploop <str>: flag for independent looping
            self.loop_delim <str>: character for delimiting loops
            self.loop_start <str>: character indicating start of loop
            self.loop_end <str>: character indicating end of loop
            self.baseinput <MASTFile>: MASTFile created from *.inp input file
            self.pegloop1 <str>: flag for one pegged loop
            self.pegloop2 <str>: flag for second pegged loop

        Looping must be indicated at the beginning of the line, and the
        text to be looped must be complete:
        indeploop mast_kpoints (3x3x3 G, 5x5x5 G, 2x2x2 M, 4x4x4 M) 
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.indeploop = "indeploop"
        self.loop_delim = ","
        self.loop_start = "("
        self.loop_end = ")"
        self.baseinput = MASTFile(self.keywords['inputfile'])
        self.pegloop1 = "pegloop1"
        self.pegloop2 = "pegloop2"
    
    def main(self, verbose=0):
        """Scan for independent loops and set up dictionaries.
            Return:
                createdfiles <list of str>: list of created input files
        """
        indepdict=self.scan_for_loop(self.indeploop)
        pegdict1 = self.scan_for_loop(self.pegloop1)
        pegdict2 = self.scan_for_loop(self.pegloop2)
        if len(indepdict.keys()) == 0 and len(pegdict1.keys()) == 0 and len(pegdict2.keys()) == 0:
            return dict()
        alldict = dict(indepdict)
        alldict.update(pegdict1)
        alldict.update(pegdict2)
        indepcomb=self.get_combo_list(indepdict, 0)
        pegcomb1=self.get_combo_list(pegdict1, 1)
        pegcomb2=self.get_combo_list(pegdict2, 1)
        allcombs = self.combine_three_combo_lists(indepcomb, pegcomb1, pegcomb2)
        datasets = self.prepare_looped_datasets(alldict, allcombs)
        createdfiles = self.create_input_files(datasets)
        if verbose == 1:
            self.print_list(indepcomb)
            self.print_list(pegcomb1)
            self.print_list(pegcomb2)
            self.print_list(allcombs)
            for datakey in datasets:
                self.print_list(datasets[datakey])
        return createdfiles

    def scan_for_loop(self, loopkey):
        """Scan for loops.
            Args:
                loopkey <str>: string for searching for loops. This is
                        either "indeploop" or "pegloop"
            Return:
                loopdict <dict>: Dictionary where loop line indices in the file
                                are keys, with values
                                'looplist', containing a list of strings, and
                                'prepend', containing any text to prepend
                                'append', containing any text to append
                    e.g. loopdict[12]['looplist']=['3x3x3 G','5x5x5 G']
                                     ['prepend']='mast_kpoints '
                                     ['append']=''
        """
        loopdict=dict()
        realline=""
        split1=""
        split2=""
        numlines = len(self.baseinput.data)
        lidx =0
        while lidx < numlines:
            dataline = self.baseinput.data[lidx].strip()
            if loopkey in dataline:
                loopdict[lidx] = dict()
                realline=dataline.split(' ',1)[1] #strip off indeploop
                split1 = realline.split(self.loop_start)
                split2 = split1[1].split(self.loop_end)
                loopdict[lidx]['prepend'] = split1[0]
                loopdict[lidx]['append'] = split2[1]
                loopdict[lidx]['looplist'] = split2[0].split(self.loop_delim)
            lidx = lidx + 1
        #print "TTM DEBUG: ", loopdict
        if len(loopdict.keys()) == 0:
            return dict()
        return loopdict
    
    def get_combo_list(self, loopdict, pegged=0):
        """Prepare a combination list of looping indices.
            Args:
                loopdict <dict>: dictionary of looping items from
                                    scan_for_loop
            Returns:
                combolist <list of str>: list of strings, like ["3-0","5-1"]
                            representing combinations
        """
        combolist=list()
        flatlists=list()
        loopkeys = list(loopdict.keys())
        loopkeys.sort()
        if pegged == 0:
            for loopkey in loopkeys:
                numloop = len(loopdict[loopkey]['looplist'])
                loopct=0
                flatlist=list()
                while loopct < numloop:
                    flatlist.append(str(loopkey) + '-' + str(loopct))
                    loopct = loopct + 1
                flatlists.append(flatlist)
            import itertools
            prod_list = itertools.product(*flatlists)
            stopiter = 0
            while not stopiter:
                try:
                    mycomb = prod_list.next()
                except StopIteration:
                    stopiter = 1
                if stopiter == 0:
                    combolist.append(list(mycomb))
        elif pegged == 1:
            if len(loopkeys) == 0:
                return combolist #Empty list
            numloop = len(loopdict[loopkeys[0]]['looplist']) #all same len
            numct=0
            while numct < numloop:
                flatlist=list()
                for loopkey in loopkeys:
                    flatlist.append(str(loopkey) + '-' + str(numct))
                numct = numct + 1
                combolist.append(flatlist)
        #print "TTM DEBUG: ", flatlists
        return combolist

    def combine_three_combo_lists(self, indeplist, peglist1, peglist2):
        """Combine two pegged lists and one independent list.
            Args:
                indeplist <list of list>: List of indeploop combinations.
                peglist1 <list of list>: List of pegged loop 1 combinations
                peglist2 <list of list>: List of pegged loop 2 combinations
            Returns:
                alllist <list of list>: List of all combinations
        """
        templist=list()
        threelist=list()
        templist = self.combine_combo_lists(indeplist, peglist1)
        threelist = self.combine_combo_lists(templist, peglist2)
        return threelist

    def combine_combo_lists(self, indeplist, peggedlist):
        """Combine combination lists.
            Args:
                indeplist <list of list>: List of indeploop combinations.
                peggedlist <list of list>: List of pegged loop combinations
            Returns:
                alllist <list of list>: List of all combinations
        """
        alllist=list()
        if len(peggedlist) == 0:
            return indeplist
        if len(indeplist) == 0:
            return peggedlist
        for pegitem in peggedlist:
            for indepitem in indeplist:
                alllistentry = list(pegitem)
                alllistentry.extend(indepitem)
                alllist.append(alllistentry)
        return alllist



    def print_list(self, mylist):
        """Print a list.
            Args:
                mylist <list of str>
        """
        for myitem in mylist:
            print myitem
        return


    def prepare_looped_lines(self, alldict, comblist):
        """Prepare looped lines from looping dictionary.
            Args:
                loopdict <dict>: dictionary of looping items from
                                    scan_for_loop
            Returns:
                loopline_dict <dict of dict>: dictionary of different lines for
                                        input files, with keys being
                                        the looped line index
                
        """
        loopline_dict=dict()
        for stridx in comblist:
            lidx = int(stridx.split('-')[0])
            loopidx = int(stridx.split('-')[1])
            loopline_dict[lidx] =  alldict[lidx]['prepend'] + alldict[lidx]['looplist'][loopidx].strip() + alldict[lidx]['append'] + '\n'
        return loopline_dict

    def prepare_looped_datasets(self, alldict, allcombs):
        """Prepare looped datasets from looping lines.
            Args:
                alldict <dict>: line dictionary for looping
                allcombs <list of list>: index combinations for looping
            Returns:
                datasets_dict <dict of list>: full datasets for new input files
        """
        datasets_dict=dict()
        numcombs = len(allcombs)
        combct = 0
        while combct < numcombs:
            newdata = list(self.baseinput.data)
            loopedlines = dict()
            loopedlines = self.prepare_looped_lines(alldict, allcombs[combct])
            for lvalidx in loopedlines.keys():
                newdata[lvalidx] = loopedlines[lvalidx]
            datasets_dict[combct] = newdata
            combct = combct + 1
        return datasets_dict

    def create_input_files(self, datasets_dict):
        """Create independently looped input files.
            Args:
                datasets_dict <dict of list>

            Returns:
                createdfiles <list of str>: list of newly-created file names
            Creates an input file for each entry.


        """
        ifname = self.keywords['inputfile']
        dirstem = os.path.dirname(ifname)
        basename = os.path.basename(ifname).split('.')[0]
        createdfiles=list()
        if dirstem == "":
            dirstem = os.getcwd()
        dkeys = datasets_dict.keys()
        dkeys.sort()
        dct=1
        for didx in dkeys:
            newfile = MASTFile()
            newfile.data = list(datasets_dict[didx])
            newname="%s/loop_%s_%s.inp" % (dirstem, basename, str(dct).zfill(2))
            newfile.to_file(newname)
            #createdfiles.append(os.path.basename(newname))
            createdfiles.append(newname)
            dct=dct+1
        return createdfiles

