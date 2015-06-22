##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import logging
from MAST.utility import loggerutils
from MAST.utility import MASTFile
from MAST.utility import MASTError

"""Utility to read the indented recipe file."""
def read_recipe(rawrlist, verbose=0):
    """Read the indented recipe.
        "Recipe" is a protected keyword signaling the recipe name.
        Args:
            rawrlist <list of str>: List of recipe lines from $recipe section of input file
        Returns:
            totpdict <dict>: Dictionary of 
                [parentname][child]=[parent method group]
            totcdict <dict>: Dictionary of
                [childname]['parents']=[parent,parent,...]
                [childname]['method']=[method group]
            rname <str>: Recipe name
    """
    logger=loggerutils.get_mast_logger("read_recipe")

    rfile = list(rawrlist) #MASTFile(filename)
    rdata = list()
    #preprocess by removing blank lines and any "recipe" line
    for line in rfile:
        duplicate = line
        if ((len(duplicate) - len(duplicate.lstrip(' '))) % 4 != 0):
            raise MASTError("recipe/recipeutility", "Recipe at %s contains incorrect number of whitespace chars at the beginning of the line! Please convert all indentations to the appropriate number of groups of four spaces." % filename)
        if '\t' in line:
            raise MASTError("recipe/recipeutility", "Recipe at %s contains tabs! Please convert all indentations to the appropriate number of groups of four spaces." % filename)
        myline = line.rstrip() #right-hand strip of carriage return only
        if len(myline) == 0: #blank line
            pass
        #elif (len(myline) > 6) and (myline[0:6].lower() == "recipe"):
            #rname = myline.split()[1]	  
        elif myline.strip()[0] == "#": #comment line, hash starting wherever
            pass
        else:
            rdata.append(myline)
    subrdict = split_into_subrecipes(rdata) 
    howtorun=dict()
    parentstocheck=dict()
    howtoupdate=dict()
    for subkey in subrdict.keys():
        idict = make_indentation_dictionary(subrdict[subkey])
        if verbose == 1:
            ikeys = idict.keys()
            ikeys.sort()
            for indentkey in ikeys:
                logger.info("Indent: %s" % indentkey)
                for myitem in idict[indentkey]:
                    logger.info(myitem)
        [onehtu, oneptc, onehtr]=parse_indentation_dict(idict)
        for pkey in onehtu.keys():
            if not (pkey in howtoupdate.keys()):
                howtoupdate[pkey]=dict()
            for pmkey in onehtu[pkey].keys():
                howtoupdate[pkey][pmkey]=onehtu[pkey][pmkey]
        for ckey in oneptc.keys():
            if not (ckey in parentstocheck.keys()):
                parentstocheck[ckey]=list()
            parentstocheck[ckey].extend(oneptc[ckey])
            howtorun[ckey]=onehtr[ckey]
    if verbose==1:
        logger.info("How-to-update-children tree: ")
        keylist = howtoupdate.keys()
        keylist.sort()
        for htukey in keylist:
            logger.info("%s: %s" % (htukey, howtoupdate[htukey]))
        clist = parentstocheck.keys()
        clist.sort()
        logger.info("Parents-to-check tree: ")
        for ckey in clist:
            logger.info("%s: %s" % (ckey, parentstocheck[ckey]))
        logger.info("How-to-run tree: ")
        for ckey in clist:
            logger.info("%s: %s" % (ckey, howtorun[ckey]))

    return [howtoupdate, parentstocheck, howtorun]

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
    """Parse an indentation dictionary into three parts. 
        Args:
            idict <dict>: Indentation dictionary
        Returns:
            how_to_update <dict>: Dictionary of how to update
                from the parent to the child:
                how_to_update['perfect_opt1']['perfect_opt2']=
                            ['volrelax_to_singlerun']
            parents_to_check <dict>: Dictionary of parents
                for the child to check (as a list)
                parents_to_check['neb_1-2_opt1']=
                        ['defect1_stat','defect2_stat']
            how_to_run <dict>: Dictionary of how to write,
                evaluate run readiness, run, and evaluate
                completion of the ingredient:
                how_to_run['perfect_opt1']=['volrelax_to_singlerun']
    """
    htu=dict()
    ptc=dict()
    htr=dict()
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
            if not (pname in htu.keys()):
                htu[pname]=dict()
            if not (iidx + 1) in ikeys: #no children
                pass
            else:
                cllist = idict[iidx+1]
                for clitem in cllist:
                    clnum = clitem[0]
                    clname = clitem[1]
                    clmethod = clitem[2]
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
                        htu[pname][clname]=pmethod
                        if not (clname in ptc.keys()):
                            ptc[clname]=list()
                            htr[clname]=clmethod 
                        ptc[clname].append(pname)
            if not (pname in ptc.keys()):
                ptc[pname]=list() #add tree originator(s)
                htr[pname]=pmethod
            pidx=pidx+1
        iidx=iidx+1
    return [htu,ptc,htr]

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

