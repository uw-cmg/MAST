import os
import sys
import fnmatch
from MAST.utility import MASTError

def walkdirs(existdir, mindepth=1, maxdepth=5, matchme=""):
    """Walk through directory and return list of subdirectories."""
    if not(os.path.exists(existdir)):
        raise MASTError("utility","No directory at " +existdir)
#   walk and make main list
    walkme=""
    walkentry=""
    walkme = os.walk(existdir)
    if walkme == ():
        raise MASTError("utility walkdirs","No folders found in " + existdir)
    bigdirlist=[]
    for walkentry in walkme:
        bigdirlist.append(walkentry[0])
#   sort list and pare down according to min/max depth
    numsep=0
    minnumsep=0
    maxnumsep=0
    numsep = existdir.count('/') # find base number of forward slashes
    minnumsep = numsep + mindepth
    maxnumsep = numsep + maxdepth
    onedir=""
    onenumsep=0
    smalldirlist=[]
    for onedir in bigdirlist:
        onenumsep = onedir.count('/')
        if (onenumsep >= minnumsep) and (onenumsep <= maxnumsep):
            smalldirlist.append(onedir)
    smalldirlist.sort()
    if not matchme == "":
        matchdirlist=[]
        for mydir in smalldirlist:
            if fnmatch.fnmatch(mydir, matchme):
                matchdirlist.append(mydir)
        return matchdirlist
    else:
        return smalldirlist

def walkfiles(existdir, mindepth=1, maxdepth=5, matchme=""):
    """Walk through directory and subdirectories and return list of files."""
    if not(os.path.exists(existdir)):
        raise MASTError("utility","No directory at " +existdir)
#   walk and make main list
    walkme=""
    walkentry=""
    walkme = os.walk(existdir)
    if walkme == ():
        raise MASTError("utility walkfiles","No folders found in " + existdir)
    filetree=""
    diritem=""
    fileitem=""
    fullfile=""
    filelist=[]
    filetree = os.walk(existdir)
    for diritem in filetree:
        for fileitem in diritem[2]:
            fullfile = diritem[0] + '/' + fileitem
            filelist.append(fullfile)
    #for fullfile in filelist:
    #    print fullfile
#   #   sort and pare file list
    paredfilelist=[]
    numsep=0
    minnumsep=0
    maxnumsep=0
    numsep = existdir.count('/') # find base number of forward slashes
    minnumsep = numsep + mindepth
    maxnumsep = numsep + maxdepth
    onefile=""
    onenumsep=0
    for onefile in filelist:
        onenumsep = onefile.count('/')
        if (onenumsep >= minnumsep) and (onenumsep <= maxnumsep):
            paredfilelist.append(onefile)
    paredfilelist.sort()
    if not matchme == "":
        matchfilelist=[]
        for myfile in paredfilelist:
            if fnmatch.fnmatch(myfile, matchme):
                matchfilelist.append(myfile)
        return matchfilelist
    else:
        return paredfilelist
