import os
try:
    from mpi4py import MPI
except:
    pass

def ofile(name):
    l=open(name,'a')
    return l

def setup_output(filename, restart, nindiv, indiv_defect_write, genealogy, 
    allenergyfile, fingerprinting, debug):
    """
        Subprogram to set up the directories and file outputs for the optimizer.
        Inputs:
            filename - String name of the run output folder
            restart - True/False boolean to determine if reading from 
                previous output folder
            nindiv - Integer number of structures to generate output
            indiv_defect_write - True/False boolean to determine if iso output 
                files should be created for Defect case
            genealogy - True/False boolean to determine whether or not 
                to create an output file for the genealogy of the run structures
            allenergies - True/False boolean to determine whether or not to 
                create an output file for all the energies of each generation
            fingerprinting - True/False boolean to determine whether or not 
                to create the output files for fingerprint distances and 
                minimum fingerprint by generation
            debug - List of strings.  If debug !=['None'], then a debug 
                file will be created
        Output:
            Dictionary containing output files where the key describes the file 
            structure. Keys include:
                Files = List of file objects for individual structures in a population
                ifiles = List of file objects for isolated sections of individual 
                    structures in a population for a Defect run
                output = File object for general output
                summary = File object for general summary
                Genealogyfile = File object for genealogy data or None
                tenergyfile = File object for all energy data or None
                debugfile = File object for debug file or None
                fpfile = File object for fingerprint distance data or None
                fpminfile = File object for minimum fingerprint data or None
            Note that with the exception of the general output file, all files will
            be stored in a folder under the filename input given
    """
    try:
        rank = MPI.COMM_WORLD.Get_rank()
    except:
        rank = 0
    path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(filename,rank))
    if not os.path.exists(path):
        os.mkdir(path)
        if restart:
            raise RuntimeError('Cannot find directory for restart', path)
    flist = [os.path.join(path,'indiv{0:02d}.xyz'.format(x)) for x in range(nindiv)]
    flist.append(os.path.join(path,'StructureSummary.txt'))
    files = map(ofile,flist)
    if indiv_defect_write:
        iflist = [os.path.join(path,'indiv{0:02d}-iso.xyz'.format(x)) for x in range(nindiv)]
        ifiles = map(ofile,iflist)
    else:
        ifiles = None
    outname =  os.path.join(os.getcwd(),'{0}-rank{1}.txt'.format(filename,rank))
    output = open(outname,'a')
    summary = open(os.path.join(path,'Summary-{0}.txt'.format(filename)),'a')
    optifile = os.path.join(path,'Optimizer-restart-file.txt')
    if genealogy:
        Genealogyfile = open(os.path.join(path,'Genealogy-{0}.txt'.format(filename)),'a')
    else:
        Genealogyfile = None
    if allenergyfile: 
        tenergyfile = open(os.path.join(path,'All_Energies-{0}.txt'.format(filename)),'a')
    else:
        tenergyfile = None
    if debug != ['None']:
        debugfile = open(os.path.join(path,'Debug.xyz'),'a')
    else:
        debugfile = None
    if fingerprinting: 
        fpfile=open(os.path.join(path,'Fingerprints-{0}.txt'.format(filename)),'a')
        fpminfile=open(os.path.join(path,'FingerprintMin-{0}.txt'.format(filename)),'a')
    else:
        fpfile = None
        fpminfile = None
    outfiles = {'files':files, 'ifiles':ifiles, 'output':output, 'summary':summary,
        'Genealogyfile':Genealogyfile, 'tenergyfile':tenergyfile, 'debugfile':debugfile,
        'fpfile':fpfile, 'fpminfile':fpminfile,'optimizerfile':optifile}
    return outfiles
