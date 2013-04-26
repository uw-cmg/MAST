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

def get_mast_install_path():
    getpath = os.getenv('MAST_INSTALL_PATH')
    if getpath == None:
        raise MASTError("utility dirutil","No path set in environment variable MAST_INSTALL_PATH")
    return getpath


def directory_is_locked(dirname):
    if os.path.isfile(dirname + "/mast.write_files.lock"):
        return True
    else:
        return False

def lock_directory(dirname):
    import time
    if directory_is_locked(dirname):
        wait_to_write(dirname)
    lockfile = open(dirname + "/mast.write_files.lock", 'wb')
    lockfile.writelines(time.ctime())
    lockfile.close()

def unlock_directory(dirname):
    if not directory_is_locked(dirname):
        raise MASTError("utility unlock_directory",
            "Tried to unlock a directory which was not locked.")
    os.remove(dirname + "/mast.write_files.lock")

def wait_to_write(dirname):
    import time
    waitmax = 1000
    waitcount = 1
    while directory_is_locked(dirname) and (waitcount < waitmax):
        time.sleep(5)
        waitcount = waitcount + 1
    if directory_is_locked(dirname):
        raise MASTError("utility wait_to_write", 
            "Timed out waiting to obtain lock on directory.")
