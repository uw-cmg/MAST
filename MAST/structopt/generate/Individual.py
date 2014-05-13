import copy
from ase import Atom, Atoms
from MAST.structopt.inp_out.write_individual import write_individual

class Individual(object):
    """Defines class object individual for use in evolution"""
    def __init__(self, data, fitness=0, index=0, history_index='0', energy=0, tenergymx=0, 
                tenergymin=0, bulki=Atoms(), bulko=Atoms(), box=Atoms(), pressure=0, volume=0, 
                force=0, purebulkenpa=0, natomsbulk=0, fingerprint=0, swaplist=[],
                vacancies=Atoms(),swaps=Atoms()):
        self.fitness=fitness
        self.index=index
        self.history_index=history_index
        self.energy=energy
        self.tenergymx=tenergymx
        self.tenergymin=tenergymin
        self.bulki=bulki
        self.bulko=bulko
        self.box=box
        self.pressure=pressure
        self.volume=volume
        self.force=force
        self.purebulkenpa=purebulkenpa
        self.natomsbulk=natomsbulk
        self.fingerprint=fingerprint
        self.swaplist=swaplist
        self.vacancies=vacancies
        self.swaps=swaps
        self.data=[data]
    
    def __getitem__(self, i):
        return self.data[i]
    
    def __setitem__(self, key, item):
        self.data[key]=item
    
    def __iter__(self):
        return self
    
    def __copy__(self):
        self.data=copy.copy(self.data)
        return self.data
    
    def __deepcopy__(self):
        self.data=copy.deepcopy(self.data)
        return self.data
    
    def __setattr__( self, attr, value ):
        super( Individual, self ).__setattr__( attr, value )
    
    def next(self):
        if self.index > len(self.data):
            raise StopIteration
        self.index=self.index+1
        return self.data[self.index]
    
    def len(self):
        return len(self.data)
    
    def duplicate(offspring):
        """Duplicates individual"""
        dupli=offspring[0].copy()
        dup=Individual(dupli)
        dup.fitness=copy.copy(offspring.fitness)
        dup.index=copy.copy(offspring.index)
        dup.history_index=copy.copy(offspring.history_index)
        dup.energy=copy.copy(offspring.energy)
        dup.tenergymx=copy.copy(offspring.tenergymx)
        dup.tenergymin=copy.copy(offspring.tenergymin)
        dup.bulki=offspring.bulki.copy()
        dup.bulko=offspring.bulko.copy()
        dup.box=offspring.box.copy()
        dup.pressure=copy.copy(offspring.pressure)
        dup.volume=copy.copy(offspring.volume)
        dup.force=copy.copy(offspring.force)
        dup.purebulkenpa=copy.copy(offspring.purebulkenpa)
        dup.natomsbulk=copy.copy(offspring.natomsbulk)
        dup.fingerprint=copy.copy(offspring.fingerprint)
        dup.swaplist=copy.deepcopy(offspring.swaplist)
        dup.vacancies=offspring.vacancies.copy()
        dup.swaps=offspring.swaps.copy()
        return dup
    
#     def write(self,indivfile):
#         write_individual(self,indivfile)
#         return
    
#     def read(self,indivfile):
#         self = read_individual(indivfile)
#         return self


