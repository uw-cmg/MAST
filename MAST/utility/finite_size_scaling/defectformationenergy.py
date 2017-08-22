#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Zhewen Song
##############################################################

import shutil, os
from collections import defaultdict
import pymatgen as pmg
from MAST.utility.defect_formation_energy.potential_alignment import PotentialAlignment  
import pandas
from MAST.utility.finite_size_scaling.parse_input import ParseInput 

class DefectFormationEnergy:
    
    def __init__(self):
        self.inputs = ParseInput().read_input()
        self.pa = PotentialAlignment()
    
    def read_CARs(self,GGA,HSE):
        perfect = defaultdict(defaultdict)
        defect = defaultdict(lambda:defaultdict(defaultdict))
        PA = dict()
        CARs = os.listdir('.')
        perfect_tag = self.inputs['perfect']
        defect_tag = self.inputs['defect']
        for i in range(len(CARs)):
            keywords = CARs[i].split('_')
            if not 'CAR' in CARs[i]: continue
            if 'HSE' in CARs[i]:
                size = keywords[1]+'_HSE'
            elif keywords[1] in ' '.join(HSE['size']):
                size = keywords[1]+'_GGA'
            else: size =  keywords[1]
            
            if keywords[0]==perfect_tag and keywords[-1]=='OSZICAR':
                perfect[size]['energy'] = float(open(CARs[i]).readlines()[-1].split('E0=')[1].split()[0])
            elif keywords[0]==defect_tag and keywords[-1]=='OSZICAR':
                defect[keywords[2]][keywords[3]][size] = float(open(CARs[i]).readlines()[-1].split('E0=')[1].split()[0])
            elif keywords[0]==perfect_tag and keywords[-1]=='OUTCAR':
                PA[size] = self.pa.read_outcar(CARs[i])
            elif keywords[0]==defect_tag and keywords[-1]=='OUTCAR':
                PA['%s_%s_%s'%(keywords[2],keywords[3],size)] = self.pa.read_outcar(CARs[i])
            elif keywords[0]==perfect_tag and (keywords[-1]=='POSCAR' or keywords[-1]=='CONTCAR'):
                shutil.copy(CARs[i],'POSCAR')
                Ele = pmg.Structure.from_file('POSCAR').species
                ele = dict()
                for k in range(len(Ele)): 
                    if not Ele[k] in ele.keys():
                        ele[Ele[k]] = 1
                    else: ele[Ele[k]] += 1
                perfect[size]['ele'] = ele
            elif keywords[0]==defect_tag and (keywords[-1]=='POSCAR' or keywords[-1]=='CONTCAR'):
                shutil.copy(CARs[i],'POSCAR')
                Ele = pmg.Structure.from_file('POSCAR').species
                ele = dict()
                os.system('rm POSCAR')
                for k in range(len(Ele)): 
                    if not Ele[k] in ele.keys():
                        ele[Ele[k]] = 1
                    else: ele[Ele[k]] += 1
                defect[keywords[2]]['ele'][size] = ele
        return [perfect,defect,PA]



    def calculate_dfe(self,GGA,HSE):
        [perfect,defect,PA] = self.read_CARs(GGA,HSE)
        Eform = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(defaultdict))))
        for N in defect.keys():
            for Q in defect[N].keys():
                if Q=='ele': continue
                if Q[-2]=='p': charge = float(Q[-1])
                elif Q[-2]=='n': charge = -float(Q[-1])
                for S in defect[N][Q].keys():
                    eshift = self.pa.get_potential_alignment(PA[S],PA['%s_%s_%s'%(N,Q,S)])
                    if not 'HSE' in S: xc = GGA
                    else: xc = HSE
                    if not '_' in S: X = 'GGA'
                    else: X = 'HSE'
                    for C in xc['chem_pot'].keys():
                        chempot = xc['chem_pot'][C]
                        Eform[X][C][N][Q][S] = defect[N][Q][S]-perfect[S]['energy']+charge*(xc['gap']['efermi']+eshift)
                        for pe in perfect[S]['ele']:
                            Eform[X][C][N][Q][S] += perfect[S]['ele'][pe]*chempot[str(pe)]
                        for de in defect[N]['ele'][S].keys():
                            Eform[X][C][N][Q][S] -= defect[N]['ele'][S][de]*chempot[str(de)]
        return Eform

    def main(self,GGA,HSE):
        Eform = self.calculate_dfe(GGA,HSE)
        for X in Eform.keys():
            for C in Eform[X].keys():
                fp = open('%s_%s.txt'%(C,X),'w')
                for N in Eform[X][C].keys():
                    fp.write(N+'\n')
                    dfe = []
                    charge = Eform[X][C][N].keys()
                    for i in range(len(charge)):
                        if charge[i][-2]=='p': charge[i]=int(charge[i][-1])
                        elif charge[i][-2]=='n': charge[i]=-int(charge[i][-1])
                    charge.sort()
                    for i in range(len(charge)):
                        if charge[i]<0: charge[i]='q=n'+str(abs(charge[i]))
                        else: charge[i]='q=p'+str(charge[i])
                    dfe = []
                    for q in range(len(charge)):
                        Q = charge[q]
                        string = ''
                        if X=='GGA':
                            if len(Eform.keys())==2:
                                size=GGA['size']+['delta_E']
                            else: size=GGA['size']
                        elif X=='HSE': size=HSE['size']
                        for S in size:
                            if S=='delta_E':
                                for i in HSE['size']:
                                    if 'HSE' in i: hse = i
                                    else: gga = i
                                try: string += '%.4f '%(Eform['HSE'][C][N][Q][hse]-Eform['HSE'][C][N][Q][gga])
                                except KeyError: string += 'NaN  '
                            else:
                                try: string += '%.4f '%(Eform[X][C][N][Q][S])
                                except KeyError: string += 'NaN  '
                        charge[q] += ':'
                        dfe.append(string.split())
                    fp.write(pandas.DataFrame(dfe,charge,size).to_string()+'\n\n')
                        
if __name__ == "__main__":
    GGA = ParseInput().read_input()['GGA']
    HSE = ParseInput().read_input()['HSE']
    DefectFormationEnergy().main(GGA,HSE)
