import pymatgen as mg
from pymatgen.core.structure import Structure

filename=raw_input('file name:')
l=raw_input('x scale:')
m=raw_input('y scale:')
n=raw_input('z scale:')

struct=mg.read_structure(filename)
Structure.make_supercell(struct,[l,m,n])
mg.write_structure(struct,"poscar_temp")

f=open("poscar_temp")
w=open("mast.inp","w")

def getinfo(line):
    line=line.strip('\n')
    data=line.split(' ')
    while 1:
        try: data.remove('')
        except: break
    return data

line=[]
elenum=[]

line.append(f.readline())
while line[len(line)-1]:
    line.append(f.readline())

i=0
while i<len(line):
    if getinfo(line[i])==[]:
        line.remove(line[i])
        i=i-1
    i=i+1

w.write('begin elementmap\n')

data=getinfo(line[5])   

for j in range(len(data)):
    w.write('X'+str(j+1)+' '+data[j]+'\n')
w.write('end\n\n')
   
w.write('begin lattice\n')

for i in range(2,5):
    data=getinfo(line[i]) 
    w.write(data[0]+' '+data[1]+' '+data[2]+'\n')

w.write('end\n\n')

data=getinfo(line[6])

for j in range(len(data)):
    elenum.append(int(data[j]))
          
w.write('begin coordinates\n')
k=8
for i in range(len(elenum)):
    for j in range(elenum[i]):
        data=getinfo(line[k])
        k=k+1
        w.write('X'+str(i+1)+' '+data[0]+' '+data[1]+' '+data[2]+'\n')
w.write('end\n')


f.close()
w.close()



