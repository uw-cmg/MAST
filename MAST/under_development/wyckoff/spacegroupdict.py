#!/usr/bin/env python
def make_spacegroup_dict():
    """Make spacegroup table from Seto Yusuke homepage:
        http://133.50.156.143/~seto/eng/
        index.php?plugin=related&page=Table%20of%20Space%20Group
        Returns:
            spdict <dict>: Dictionary of the format
            ['hall_to_serial'][hall symbol <str>] = serial number <int>
            ['serial_to_hall'][serial number <int>] = hall symbol <str>
    """
    from MAST.utility import MASTFile

    spfile=MASTFile("spacegrouptable.txt")
    mydata=list(spfile.data)
    mydata.pop(0) #remove header line
    spdict=dict()
    spdict['hall_to_serial']=dict()
    spdict['serial_to_hall']=dict()
    myhall=""
    myser=0
    lct=0
    for line in mydata:
        lct=lct+1
        lsplit = line.strip().split("|")
        try:
            myhall = lsplit[4].strip()
            myser = int(lsplit[0].strip())
        except IndexError:
            print "Error in line %1i" % lct
            continue
        if myhall in spdict['hall_to_serial'].keys():
            print "Duplicate hall key %s for line %1i" % (myhall, lct)
            print "Existing entry is serial number %1i" % spdict['hall_to_serial'][myhall]
        else:
            spdict['hall_to_serial'][myhall]=myser
        if myser in spdict['serial_to_hall'].keys():
            print "Duplicate serial key %1i for line %1i" % (myser, lct)
            print "Existing entry is hall symbol %s" % spdict['serial_to_hall'][myser]
        else:
            spdict['serial_to_hall'][myser]=myhall
    return spdict

def make_wyckoff_dict():
    """Make wyckoff table from Seto Yusuke homepage:
        http://133.50.156.143/~seto/eng/
        index.php?plugin=related&page=Table%20of%20Space%20Group
        Returns:
            wkdict <dict>: Dictionary of the format
            [serial # <int>]['spacegroup'] = space group <str>
                            ['letters'][letter <str>]['multiplicity']= <int>
                            ['letters'][letter <str>]['positions']=<list of str>
            Note that the number of positions does not necessarily equal
            the multiplicity, e.g. (0,0,z) with multiplicity of 3 also
            includes, but does not list, (x,0,0) and (0,y,0)
    """
    from MAST.utility import MASTFile

    wkfile=MASTFile("wyckofftable.txt")
    mydata=list(wkfile.data)
    mydata.pop(0) #remove header line
    wkdict=dict()
    ltype=""
    myspace=""
    myserial=0
    mymult=0
    myletter=""
    lct=0
    for line in mydata:
        lct=lct+1
        lsplit = line.strip().split("|")
        if (lsplit[0].strip() == "") and not (lsplit[2].strip() == ""):
            ltype="letter"
        elif (lsplit[0].strip() == "") and (lsplit[2].strip() == ""):
            ltype="extrapoints"
        else:
            ltype="header"
        if ltype == "header":
            myspace = lsplit[0].strip()
            myserial = int(lsplit[1].strip())
            myletter=""
            mymult=0
            if myserial in wkdict.keys():
                print "Duplicate serial number %1i for line %1i." % (myserial, lct)
            else:
                wkdict[myserial]=dict()
                wkdict[myserial]['spacegroup']=myspace
                wkdict[myserial]['letters']=dict()
        elif ltype == "letter":
            mymult = int(lsplit[2].strip())
            myletter = lsplit[3].strip()
            if myletter in wkdict[myserial]['letters'].keys():
                print "Duplicate letter key %s for line %1i, serial number %1i. Appending to existing list." % (myletter, lct, myserial)
                for idx in range(5, len(lsplit)):
                    trypt = lsplit[idx].strip()
                    if not trypt == "":
                        wkdict[myserial]['letters'][myletter]['positions'].append(trypt)
            else:
                wkdict[myserial]['letters'][myletter]=dict()
                wkdict[myserial]['letters'][myletter]['multiplicity']=mymult
                wkdict[myserial]['letters'][myletter]['positions']=list()
                for idx in range(5, len(lsplit)):
                    trypt = lsplit[idx].strip()
                    if not trypt == "":
                        wkdict[myserial]['letters'][myletter]['positions'].append(trypt)
        elif ltype == "extrapoints":
            for idx in range(5, len(lsplit)):
                trypt = lsplit[idx].strip()
                if not trypt == "":
                    wkdict[myserial]['letters'][myletter]['positions'].append(trypt)
    return wkdict

def find_letters_for_max_multiplicity(wkdict, serialno, maxmult=10):
    """Find the number of wyckoff letters that equals or just exceeds
        a certain multiplicity threshold:
        e.g sites 4a, 2b, 2c, 4d, 2e,... with maxmult=10 would return
        letters 'a','b','c','d', for a total multiplicity of 12.
        Args:
            wkdict <dict>: Wyckoff position dictionary, from make_wyckoff_dict
            serialno <int>: serial number
            maxmult <int>: Maximum multiplicity (default 10)
        Returns:
            letterlist <list of str>: List of letters.
    """
    letterkeys = wkdict[serialno]['letters'].keys()
    letterkeys.sort()
    totmult=0
    mymult=0
    letterlist=list()
    for letter in letterkeys:
        mymult = wkdict[serialno]['letters'][letter]['multiplicity']
        totmult = totmult + mymult
        letterlist.append(letter)
        if totmult >= maxmult:
            break
    return letterlist

myspdict=make_spacegroup_dict()
mywkdict=make_wyckoff_dict()
print mywkdict[523]
myletterlist=find_letters_for_max_multiplicity(mywkdict,523,20)
print myletterlist


