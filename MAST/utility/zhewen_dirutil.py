##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys
import fnmatch
import time
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata

def immediate_subdirs(existdir, verbose=0):
    """Walk through directory and return immediate subdirectories
        (only one level)
        as directory names (not full paths)
        Args:
            existdir <str>: Directory, which should exist
            verbose <int>: 0 - concise (default)
                           1 - verbose
    """
    if not (os.path.exists(existdir)):
        raise MASTError("utility","No directory at " +existdir)
    listdir = os.listdir(existdir)
    subdirs=list()
    for myentry in listdir:
        trydir = os.path.join(existdir, myentry)
        if os.path.isdir(trydir):
            subdirs.append(myentry)
    subdirs.sort()
    if verbose == 1:
        print "subdirectories:", subdirs
    return subdirs


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
    """Walk through directory and subdirectories and return list of files.
        Args:
            existdir <str>: directory under which to search.
            mindepth <int>: minimum depth of folders to search;
                            default 1 = within that directory
            maxdepth <int>: maximum depth of folders to search; default 5
            matchme <str>: string to match; every file ending in matchme
                            will be found, since a * is required at the front
                            to match the full paths.
        Returns:
            <list of str>: list of full file paths
    """
    if not(os.path.exists(existdir)):
        raise MASTError("utility","No directory at " +existdir)
#   walk and make main list
    walkme=""
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
        if not matchme[0] == "*":
            matchme="*"+matchme
        for myfile in paredfilelist:
            if fnmatch.fnmatch(myfile, matchme):
                matchfilelist.append(myfile)
        return matchfilelist
    else:
        return paredfilelist

def get_mast_install_path():
    import MAST
    return os.path.dirname(MAST.__file__)
def get_mast_scratch_path():
    getpath = os.getenv('MAST_SCRATCH')
    if getpath == None:
        raise MASTError("utility dirutil","No path set in environment variable MAST_SCRATCH")
    return getpath

def get_mast_archive_path():
    getpath = os.getenv('MAST_ARCHIVE')
    if getpath == None:
        raise MASTError("utility dirutil","No path set in environment variable MAST_ARCHIVE")
    return getpath

def get_mast_control_path():
    getpath = os.getenv('MAST_CONTROL')
    if getpath == None:
        raise MASTError("utility dirutil","No path set in environment variable MAST_CONTROL")
    return getpath
def get_mast_recipe_path():
    raise NotImplementedError("MAST_RECIPE_PATH is an obsolete environment variable and should no longer be used.")
    getpath = os.getenv('MAST_RECIPE_PATH')
    if getpath == None:
        raise MASTError("utility dirutil","No path set in environment variable MAST_RECIPE_PATH")
    return getpath

def get_mast_platform():
    getplatform = os.getenv('MAST_PLATFORM')
    if getplatform == None:
        raise MASTError("utility dirutil","No platform set in environment variable MAST_PLATFORM")
    return getplatform.strip().lower()

def directory_is_locked(dirname):
    if os.path.isfile(dirname + "/mast.write_files.lock"):
        return True
    else:
        return False

def lock_directory(dirname, waitmax=10):
    """Lock a directory using a lockfile.
        Args:
            dirname <str>: Directory name
            waitmax <int>: maximum number of 5-second waits
    """
    if directory_is_locked(dirname):
        wait_to_write(dirname, waitmax)
    lockfile = open(dirname + "/mast.write_files.lock", 'wb')
    lockfile.writelines(time.ctime())
    lockfile.close()

def unlock_directory(dirname, waitmax=10):
    """Unlock a directory by removing the lockfile.
        Args:
            dirname <str>: Directory name
            waitmax <int>: maximum number of 5-second waits
    """
    if directory_is_locked(dirname):
        os.remove(dirname + "/mast.write_files.lock")
    else:
        waitct=1
        okay=0
        while (not directory_is_locked(dirname)) and (waitct <= waitmax):
            if directory_is_locked(dirname):
                os.remove(dirname + "/mast.write_files.lock")
                okay=1
                break
            time.sleep(5)
            waitct=waitct+1
        if not okay==1:
            raise MASTError("utility unlock_directory",
                "Tried to unlock directory %s which was not locked." % dirname)

def wait_to_write(dirname, waitmax=10):
    """Wait to write to directory.
        Args:
            dirname <str>: Directory name
            waitmax <int>: maximum number of 5-second waits
    """
    if waitmax < 1:
        waitmax = 1
    waitcount = 1
    while directory_is_locked(dirname) and (waitcount < waitmax):
        time.sleep(5)
        waitcount = waitcount + 1
    if directory_is_locked(dirname):
        raise MASTError("utility wait_to_write", 
            "Timed out waiting to obtain lock on directory %s" % dirname)
def search_for_metadata_file(metastring="",dirname="", metafilename="metadata.txt", verbose=0):
    """Match a metadata file based on input.
        Args:
            metastring <str>: equals-sign-separated metatag=value pairs with
                                commas separating the meta sections.
                Example: "ingredtype=phonon, neblabel=vac1-vac2, charge=0"
            dirname <str>: directory name to start. Default "" goes to ARCHIVE.
            metafilename <str>: metadata file name. Default "metadata.txt"
            verbose <int>: 1 for verbose messages, 0 otherwise
        Returns:
            dlist <list of str>: list of directories containing matching
                                    metadata files.
    """
    if dirname=="":
        dirname = get_mast_archive_path()
    allmetas = walkfiles(dirname, 1, 5, metafilename)
    if len(allmetas) == 0:
        raise MASTError("utility dirutil, search_for_metadata_file", "No matching metafiles found in %s for tags %s." % (dirname, metastring))
    metaparse=dict()
    metasplit = metastring.split(",")
    for metaitem in metasplit:
        onesplit=metaitem.strip().split("=")
        metaparse[onesplit[0].strip()]=onesplit[1].strip()
    metamatch=list()
    mustmatch=len(metaparse.keys())
    for mtry in allmetas:
        mokay=0
        mymeta = Metadata(metafile=mtry)
        for metatag,metaval in metaparse.iteritems():
            searchresult = mymeta.search_data(metatag)
            if searchresult[1] == metaval:
                mokay=mokay+1
            else:
                if verbose == 1:
                    print metatag, searchresult
        if mokay == mustmatch:
            metamatch.append(mtry)
    if verbose==1:
        print allmetas
        print metaparse
        print metamatch
    return metamatch

def list_methods(myclass=None, printout=1):
    """List the methods in a class.
        Args:
            myclass <Class>: Class, like BaseIngredient
            printout <int>: Print the list to screen
                            1 - print (default)
                            0 - do not print to screen
    """
    import inspect
    mylist=inspect.getmembers(myclass, predicate=inspect.ismethod)
    parsedlist=list()
    for myentry in mylist:
        parsedlist.append(myentry[0])
    if printout > 0:
        for myentry in parsedlist:
            print myentry
    return parsedlist

def dir_is_in_scratch(mydir):
    """Check if the directory is in the $MAST_SCRATCH directory.
        Args:
            mydir <str>: directory name (not full path)
        Returns:
            True or False
    """
    scratch = get_mast_scratch_path()
    dirs_in_scratch = os.listdir(scratch)
    if (mydir in dirs_in_scratch):
        return True
    else:
        return False
def get_test_dir(testname):
    """Get testing directory name.
        Args:
            testname <str>: test directory name, e.g. test_mod_1
        Returns:
            <str>: full test directory, depending on where python is
                    installed and where the tests are located, e.g.
        <python install dir>/lib/python2.7/site-packages/MAST/test/test_mod_1
        where <python install dir> might be:
                //home/username/.local/
                //share/apps/EPD_64bit
                or some other directory
        Assumes that the user starts running the nosetests in the MAST
        directory or in MAST/test
    """
    testdir=""
    curdir = os.getcwd()
    curdir = os.path.normpath(curdir)
    if os.path.basename(curdir) == testname:
        return curdir
    testdir = curdir
    if os.path.basename(curdir) == "MAST":
        testdir = os.path.join(testdir, "test")
    testdir = os.path.join(testdir, testname)
    return testdir

def get_version():
    """Get MAST version from MAST/_version.py file"""
    vfile = open(os.path.join(get_mast_install_path(),"_version.py"),"rb")
    vline = vfile.readline()
    vnum = vline.split("=")[1].strip().strip('"')
    return vnum
