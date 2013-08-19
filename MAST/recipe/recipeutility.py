############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os
from MAST.utility import MASTFile

"""Utility to read the indented recipe file."""
def read_recipe(filename, verbose=1):
    """Read the indented recipe.
        "Recipe" is a protected keyword signaling the recipe name.
        Args:
            filename <str>: Full path to recipe file
        Returns:
            totpdict <dict>: Dictionary of 
                [parentname][parent method group]=[child,child..]
            totcdict <dict>: Dictionary of
                [childname]=[parent,parent,...]
            rname <str>: Recipe name
    """
    rfile = MASTFile(filename)
    rdata = list()
    #preprocess by removing blank lines and any "recipe" line
    for line in rfile.data:
        myline = line.rstrip() #right-hand strip of carriage return only
        if len(myline) == 0: #blank line
            pass
        elif (len(myline) > 6) and (myline[0:6].lower() == "recipe"):
            rname = myline.split()[1]
        elif myline.strip()[0] == "#": #comment line, hash starting wherever
            pass
        else:
            rdata.append(myline)
    subrdict = split_into_subrecipes(rdata) 
    totpdict=dict()
    totcdict=dict()
    for subkey in subrdict.keys():
        idict = make_indentation_dictionary(subrdict[subkey])
        if verbose == 1:
            ikeys = idict.keys()
            ikeys.sort()
            for indentkey in ikeys:
                print "Indent: ", indentkey
                for myitem in idict[indentkey]:
                    print myitem
        [onepdict,onecdict]=parse_indentation_dict(idict)
        for pkey in onepdict.keys():
            if not (pkey in totpdict.keys()):
                totpdict[pkey]=dict()
            for pmkey in onepdict[pkey].keys():
                if not (pmkey in totpdict[pkey].keys()):
                    totpdict[pkey][pmkey]=list()
                totpdict[pkey][pmkey].extend(onepdict[pkey][pmkey])
        for ckey in onecdict.keys():
            if not (ckey in totcdict.keys()):
                totcdict[ckey]=list()
            totcdict[ckey].extend(onecdict[ckey])
    plist = totpdict.keys()
    plist.sort()
    if verbose==1:
        for pkey in plist:
            print pkey, ":", totpdict[pkey]
        clist = totcdict.keys()
        clist.sort()
        for ckey in clist:
            print ckey, ":", totcdict[ckey]
    return [totpdict, totcdict, rname]

def split_into_subrecipes(mydata):
    """Split an entire recipe into subrecipes. Each new zero-level 
        indentation starts a new subrecipe.
        Args:
            mydata <list>: List of lines
        Returns:
            subrecipes <dict>: Dictionary of subrecipe lines (lists)
                subrecipe[1]=list()
    """
    subrecipes=dict()
    ridx=0
    for myline in mydata:
        if get_indent_level(myline) == 0:
            ridx = ridx + 1
            subrecipes[ridx]=list()
        subrecipes[ridx].append(myline)
    return subrecipes

def get_indent_level(myline):
    """Get the indentation level of a line.
        Every 4 spaces is an indentation level.
        Args:
            line <str>: Line to get the indentation level of.
        Returns:
            ilevel <int>: Indentation level. 0 = "writing"
                                             1 = "    writing"
                                             2 = "        writing" 
                                             etc.
                            -1 if indentation level was not determined.
    """
    idx=0
    ilevel=0
    mylen = len(myline)
    while idx < mylen:
        if not (myline[idx] == " "):
            if myline[0:idx].strip() == "": #all previous chars are spaces
                return ilevel
        idx = idx + 4
        ilevel = ilevel + 1
    return -1

def parse_for_name_and_instructions(myline):
    """Parse line for name(s) and instruction name(s).
        Args: 
            myline <str>: Line to parse. Should look something like:
            defect1_stat (singlerun_to_phonon)
            OR
            defect1_stat (singlerun_to_neb), defect2_stat (singlerun_to_neb)
        Returns:
            pdict <dict>: Dictionary of name = instructions
    """
    mysplit = myline.split(',')
    pdict=dict()
    for ppiece in mysplit:
        ppiecestr = ppiece.strip()
        ppsplit = ppiecestr.split()
        pname = ppsplit[0]
        if len(ppsplit) > 1:
            pmethod = ppsplit[1][1:-1] #remove parentheses
        else:
            pmethod = "default"
        pdict[pname] = pmethod
    return pdict

def parse_indentation_dict(idict):
    """Parse an indentation dictionary 
        Args:
            idict <dict>: Indentation dictionary
        Returns:
            parentdict <dict>: Dictionary of parents and 
                instructions, with values as list of applicable 
                children for those sets of instructions
                    parentdict['perfect_opt1']
                    ['volrelax_to_singlerun']=['perfect_opt2']
            childdict <dict>: Dictionary of children, 
                        with values as list of parents
    """
    parentdict=dict()
    childdict=dict()
    iidx=0
    ikeys = idict.keys()
    ikeys.sort()
    while iidx in ikeys:
        plen = len(idict[iidx])
        pidx = 0
        while pidx < plen:
            lct = idict[iidx][pidx][0]
            pname = idict[iidx][pidx][1]
            pmethod = idict[iidx][pidx][2]
            if not (pname in parentdict.keys()):
                parentdict[pname]=dict()
            if not (pmethod in parentdict[pname].keys()):
                parentdict[pname][pmethod]=list()
            if not (iidx + 1) in ikeys: #no children
                pass
            else:
                cllist = idict[iidx+1]
                for clitem in cllist:
                    clnum = clitem[0]
                    clname = clitem[1]
                    addme=0
                    if (clnum > lct):
                        if (pidx + 1) >= plen: #no more parents
                            addme=1
                        else:
                            if lct == idict[iidx][pidx+1][0]: #dual parent
                                addme=1
                            elif clnum < idict[iidx][pidx+1][0]:
                                addme=1
                    if addme == 1:
                        parentdict[pname][pmethod].append(clname)
                        if not (clname in childdict.keys()):
                            childdict[clname]=list()
                        childdict[clname].append(pname)
            if not (pname in childdict.keys()):
                childdict[pname]=list() #add tree originator(s)
            pidx=pidx+1
        iidx=iidx+1
    return [parentdict,childdict]

def make_indentation_dictionary(subrlist):
    """Parse a single subrecipe into an indentation dictionary:
        idict[indentlevel]=list of [linenumber,name,instructions]
        Args:
            subrlist <list>: List of subrecipe lines
        Returns:
            idict <dict>: Dictionary of indentations
    """
    idict=dict() #creates name/instruction dictionaries 
                        #by indentation level
    lct=0
    for line in subrlist:
        lct=lct+1
        indent = get_indent_level(line)
        if not (indent in idict.keys()):
            idict[indent]=list()
        pdict = parse_for_name_and_instructions(line)
        for parent in pdict.keys():
            idict[indent].append([lct, parent, pdict[parent]])
    return idict

