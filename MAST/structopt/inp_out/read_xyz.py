try:
    from ase import Atom, Atoms
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt read_xyz.py, ASE must be installed."
def read_xyz(fileobj,n=-1,data=False):
    """
    Function to read multi-atom xyz file with data
    Inputs:
        fileobj = String containing file name or file object
        n = Integer indicating number for structure from xyz file
            to read. Default is last structure in file.
            Or String=='All' which will output all structures in file
        data = True/False boolean indicating whether or not to output 
            any data strings read from xyz file
    Outputs:
        atmslist = ASE Atoms class object containing structure from file
            or list of Atoms objects if n=='All'
        datalist(Optional) = String of data read from xyz file or list 
            of data strings read from xyz file
    """
    if isinstance(fileobj, str):
        fileobj = open(fileobj,'r')
    lines = fileobj.readlines()
    atmslist = []
    datalist = []
    while len(lines)>0:
        natoms = int(lines[0])
        datalist.append(lines[1])
        atm1 = Atoms()
        i =- 1
        for i in range(natoms):
            a = lines[i+2].split()
            sym = a[0]
            position = [float(a[1]),float(a[2]),float(a[3])]
            atm1.append(Atom(symbol=sym,position=position))
        atmslist.append(atm1)
        lines = lines[i+3::]
    if n == 'All':
        if data == True:
            return atmslist, datalist
        else:
            return atmslist
    else:
        if data == True:
            return atmslist[n],datalist[n]
        else:
            return atmslist[n]

