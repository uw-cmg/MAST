try:
    from ase import Atom, Atoms
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt read_xyz.py, ASE must be installed."

def read_xyz(fileobj,n=-1,data=False):
    if isinstance(fileobj, str):
        fileobj = open(fileobj,'r')
    lines = fileobj.readlines()
    fileobj.close()
    natmslist = []
    atmslist = []
    datalist = []
    readnatoms = 0
    while readnatoms<len(lines):
        natoms = int(lines[readnatoms])
        natmslist.append(natoms)
        readnatoms += natoms+2
    if n != 'All':
       if n < 0:
          n = len(natmslist)+n
       lines = lines[sum(natmslist[0:n])+2*n::]
       atm1 = Atoms()
       size = lines[1].split()
       if len(size) > 1:
          atm1.set_cell([size[1],size[2],size[3]])
       for i in range(natmslist[n]):
          a = lines[i+2].split()
          sym = a[0]
          position = [float(a[1]),float(a[2]),float(a[3])]
          atm1.append(Atom(symbol=sym,position=position))
       if data == True:
            return atm1,lines[1]
       else:
            return atm1
       return atm1
    elif n=='All':
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
      if data == True:
            return atmslist, datalist
      else:
            return atmslist


