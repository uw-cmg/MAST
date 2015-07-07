import numpy as np

class ParseInput:
    def __init__(self,inputfile=''):
        self.inputfile = inputfile
        self.method = 'L1L3'
        self.perfect = 'perfect'
        self.defect = 'defect'
        self.HSE = {'tag':'HSE','size':[],'chem_pot':dict()}
        self.GGA = {'tag':'GGA','size':[],'chem_pot':dict()}
    
    def read_gap_from_doscar(self,doscar):
        from pymatgen.electronic_structure.dos import Dos
        fp = open(doscar)
        lines = fp.readlines()
        energies=[]
        dos=[]
        efermi = float(lines[5].split()[3])
        if len(lines[6].split())==5: # spin up and dn
            spin = 2
        elif len(lines[6].split())==3: spin = 0 # no spin
        for i in range(6,len(lines)):
            energies.append(float(lines[i].split()[0]))
            if spin==2:
                dos.append(float(lines[i].split()[1])+float(lines[i].split()[2]))
            else: dos.append(float(lines[i].split()[1]))
        return {'efermi':efermi,'mygap':np.array([Dos(efermi,energies,{'spin':dos}).get_gap(spin='spin')])}
    
    def read_gap_from_eigenval(self,outcar):
        fp = open(outcar)
        lines = fp.readlines()
        cband = []
        vband = []
        for i in range(len(lines)):
            if 'band No.  band energies     occupation' in lines[i]:
                j = i+1
                while not lines[j]=='\n':
                    if float(lines[j].split()[2])==0.0 and not float(lines[j-1].split()[2])==0.0:
                        cband.append(float(lines[j].split()[1]))
                        vband.append(float(lines[j-1].split()[1]))
                        break
                    j += 1
            if 'E-fermi' in lines[i]:
                efermi = float(lines[i].split()[2])
        return {'efermi':efermi,'mygap':np.array([min(cband)-max(vband)])}
    def read_input(self):
        fp = open('dfe_input.txt', 'r')
        lines = fp.readlines()
        for l in range(len(lines)):
            keywords = lines[l].split()[0]
            if 'method' in keywords.lower():
                self.method = lines[l].split()[1]
            elif 'perfect' in keywords.lower():
                self.perfect = lines[l].split()[1]
            elif 'defect' in keywords.lower():
                self.defect = lines[l].split()[1]
            elif 'HSE' in keywords.upper():
                if 'gap' in keywords.lower():
                    if 'dos' in lines[l].split()[1].lower():
                        doscar = '%s_DOSCAR'%lines[l].split()[1]
                        GAP = self.read_gap_from_doscar(doscar)
                    elif 'band' in lines[l].split()[1].lower():
                        outcar = '%s_OUTCAR'%lines[l].split()[1]
                        GAP = self.read_gap_from_eigenval(outcar)
                    self.HSE['gap'] = GAP
                    self.HSE['gap']['mygap'] = np.append(self.HSE['gap']['mygap'], np.array(lines[l].split()[2:],dtype=float))
                elif 'size' in keywords.lower():
                    self.HSE['size'] = lines[l].split()[1:]
                else:
                    chem_pot = lines[l].split()[1:]
                    C = keywords.split('_')[0]
                    self.HSE['chem_pot'][C] = dict()
                    for ele in range((len(chem_pot))/2):
                        self.HSE['chem_pot'][C][chem_pot[2*ele]] = float(chem_pot[2*ele+1])
            else:
                if 'gap' in keywords.lower():
                    if 'dos' in lines[l].split()[1].lower():
                        doscar = '%s_DOSCAR'%lines[l].split()[1]
                        GAP = self.read_gap_from_doscar(doscar)
                    elif 'band' in lines[l].split()[1].lower():
                        outcar = '%s_OUTCAR'%lines[l].split()[1]
                        GAP = self.read_gap_from_eigenval(outcar)
                    self.GGA['gap'] = GAP
                    self.GGA['gap']['mygap'] = np.append(self.GGA['gap']['mygap'], np.array(lines[l].split()[2:],dtype=float))
                elif 'size' in keywords.lower():
                    self.GGA['size'] = lines[l].split()[1:]
                else:
                    chem_pot = lines[l].split()[1:]
                    C = keywords.split('_')[0]
                    self.GGA['chem_pot'][C] = dict()
                    for ele in range((len(chem_pot))/2):
                        self.GGA['chem_pot'][C][chem_pot[2*ele]] = float(chem_pot[2*ele+1])

        return {'method':self.method,'perfect':self.perfect,'defect':self.defect,'HSE':self.HSE,'GGA':self.GGA}
