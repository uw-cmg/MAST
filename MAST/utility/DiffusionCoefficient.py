#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Zhewen Song
# Last updated: 2015-05-xx
# FCC dilute Five-frequency model equation from R. E. Howard and J. R. Manning, Physical Review, Vol. 154, 1967.
# HCP dilute Eight-frequency model equation from P. B. Ghate, Physical Review, Vol. 133, 1963.
# BCC dilute Nine-frequency model equation from A. D. Le Claire, Physical Chemistry: An Advanced Treatise, (Academic Press, 1970)
# FCC concentrated 14-frequency model equation from Bocquet J.-L. and Le Claire A. D.
#   Bocquet (X1 = X2 = f0)
#   b1 = -18 + (w4/w3)*[(4*w1 + 14*w3)/w0]
#   b2 = (w23/w24)*(w4/w3)*[(w4/w3)^2*(8*w13+2*w11+10*w33-16*w1-28*w3)/w0 + (w4/w3)*(8*(w14-w1)+2*(w34-w3))/w0 +(10*w44-37*w4)/w0 + 63] + (w4/w3)^2*(4*w11+20*w13+11*w33-28*w1-49*w3)/w0 + (w4/w3)*(20*(w14-w1)+47*(w34-w3))/w0 + 11*w44/w0 - 173*w4/w0 + 192
#   B1 = -6 - 12*[(w23/w24)*(w4/w3)] + (w4/w3)*[4*(w21/w2)*(w23/w24)*(w4/w3) + 14*(w23/w2)]
#   D_host(c) = D_host(0) * [ 1 + b1*c + b2*c^2 + ...]
#   D_solute(c) = D_solute(0) * [ 1 + B1*c + B2*c^2 + ...]
#   c = solute concentration [atomic fraction]
##############################################################

import sys, getopt, os
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from MAST.utility import fileutil
  
class ParsingInputFiles(object):
    """Parsing files of energy and frequency.
        Attributes:
            self.inp <str>: input file name
            self.get_item_name <dict of dict>: get directories from input
            self.get_lattice <dict>: lattice information
            self.get_HB_and_HVf <float>: vacancy formation energy or/and binding energy
            self.get_barrier <dict>: energy barrier
            self.get_v <dict>: attempt frequency 
    """
    def __init__(self,inp):
        self.inp = inp
        
    def getinfo(self,line):
        """Splitting a certain line with space and enter.
            Args:
                line <str>: a certain line to be handled
        """
        line=line.strip('\n')
        data=line.split(' ')
        while 1:
            try: data.remove('')
            except ValueError: break
        return data
        
    def get_item_name(self,keyword):
        """Obtaining directories from different items such as E* and H* (the energy), v* (the frequency), frequency type (5 or 8) and lattice parameters.
            Args:
                keyword <str>: keywords that specify the desired items, including 'E', 'H', 'v', 'type', 'lattice'
        """
        item_name = dict()
        with open(self.inp,'r') as contentfile:
            content = contentfile.readlines()
        Elist=['E0','E1','E2','E3','E4','E3p','E4p','E3pp','E4pp','E5','E6','Ea','Eb','Ec','EX','Eap','Ebp','Ecp','EXp','E11','E12','E13','E14','E23','E24','E33','E34','E44']
        vlist=['v0','v1','v2','v3','v4','v3p','v4p','v3pp','v4pp','v5','v6','va','vb','vc','vX','vap','vbp','vcp','vXp','v11','v12','v13','v14','v23','v24','v33','v34','v44']
        Hlist=['HB','HVf']
        for i in range(len(content)):
            line = self.getinfo(content[i])
            if not len(line)==0:
                if (keyword=='E' and line[0] in Elist) or (keyword=='v' and line[0] in vlist):
                    item_name[line[0]] = []
                    if len(line)==2:
                        try:
                            item_name[line[0]].append(float(line[1]))
                        except ValueError:
                            pass
                    else:
                        for index in range(1,len(line)):
                            item_name[line[0]].append(line[index])
                elif (keyword=='H' and line[0] in Hlist):
                    if line[0] == 'HB':
                        if len(line) > 2:
                            item_name['HB'] = dict()
                            item_name['HB']['perfect'] = line[1]
                            item_name['HB']['vac-sub'] = line[2]
                            item_name['HB']['sub'] = line[3]
                            item_name['HB']['vac'] = line[4]
                        else:
                            item_name['HB'] = float(line[1])
                    elif line[0] == 'HVf':
                        if len(line) > 2:
                            item_name['HVf'] = dict()
                            item_name['HVf']['perfect'] = line[1]
                            item_name['HVf']['vac'] = line[2]
                        else:
                            item_name['HVf'] = float(line[1])
                elif keyword=='temp' and 'temp' in line[0]:
                    item_name['temp'] = [float(line[1]),float(line[2]),float(line[3])]
                elif keyword=='type' and 'type' in line[0]:
                    if line[0]=='type' and line[1]=='fcc_5freq':
                        item_name = 5
                    elif line[0]=='type' and line[1]=='hcp_8freq':
                        item_name = 8
                    elif line[0]=='type' and line[1]=='bcc_9freq':
                        item_name = 9
                    elif line[0]=='type' and line[1]=='fcc_14freq':
                        item_name = 14
                elif keyword=='lattice' and 'lattice' in line[0]:
                    item_name = line[1]
                elif keyword=='plotdisplay' and 'plotdisplay' in line[0]:
                    item_name = line[1]
        return item_name
    
    def get_latt(self,model,latt):
        """Obtaining lattice information from the POSCAR.
            Args:
                types <dict>: type of frequency model, types={'types':<int>}
                lattice <str>: specific name of the POSCAR given by the user
        """
        data={'a':0,'c':0,'No.':0}
        os.system('cp '+latt+'_POSCAR POSCAR')

        struct = Poscar.from_file('POSCAR').structure
        os.system('rm POSCAR')
        reduced = SpacegroupAnalyzer(struct,0.001).get_primitive_standard_structure()   
        data['No.'] = len(struct)
        if model==5:
            data['a'] = reduced.lattice.abc[0]*np.sqrt(2)*10**(-8)
        if model==9: #for BCC, primitive param is half of cubic body diagonal
            data['a'] = reduced.lattice.abc[0]*10**(-8)*2.0/np.sqrt(3.0)
        if model==14:
            data['a'] = reduced.lattice.abc[0]*np.sqrt(2)*10**(-8)
        if model==8:
            data['a'] = reduced.lattice.abc[0]*10**(-8)
            data['c'] = reduced.lattice.abc[2]*10**(-8)
        return data
        
    def get_HB_and_HVf(self,Hdir,numatom,keyword):
        """Obtaining the value of binding energy HB and vacancy formation energy HVf.
            Args:
                Hdir <dict of dict>: Hdir={'HB':{'perfect': *},{'vac-sub': *},{'sub': *},{'vac': *}},'HVf':{{'perfect': *},{'vac': *}}}, * is <file_name> or <value>
                numatom <int>: number of atoms in the system
                keyword <str>: Either 'HB' or 'HVf'
        """
        try: 
            float(Hdir[keyword])
            return Hdir[keyword]
        except TypeError:
            ene = {}
            for key in Hdir[keyword].keys():
                try: 
                    ene[key]=float(Hdir[keyword][key])
                except ValueError: 
                    ene[key] = Hdir[keyword][key]
                    with open(ene[key]+'_OSZICAR','r') as efile:
                        line = efile.readlines()
                    pt = -1   
                    while not 'E0' in line[pt]: pt = pt - 1
                    ene[key] = float(self.getinfo(line[pt])[4])
            if keyword=='HVf':
                return ene['vac'] - (numatom - 1)*ene['perfect']/numatom
            elif keyword=='HB': 
                return ene['perfect'] + ene['vac-sub'] - ene['sub'] - ene['vac']
    
    def get_barrier(self,Edir,Edir_saddle,Edir_min):
        """Obtaining the energy barriers.
            Args:
                Edir <list of dict>: Edir={'E0':*,'E1':*,...}, * is the list like [<local_min>,<saddle>]
                Edir_saddle <dict>: Edir_saddle={'E0':<saddle>,...}
                Edir_min <dict>: Edir_min={'E0':<local_min>,...}
        """        
        enebarr = dict()
        enesaddle = dict()
        eneend = dict()
        for freq in Edir.keys():
            if not len(Edir[freq])==1: 
                try: enesaddle[freq]=float(Edir_saddle[freq])
                except ValueError:
                    with open(Edir_saddle[freq]+'_OSZICAR','r') as esfile:
                        line = esfile.readlines()
                    pt = -1
                    while not 'E0' in line[pt]: pt = pt - 1
                    enesaddle[freq] = float(self.getinfo(line[pt])[4])
        for freq in Edir.keys():
            if not len(Edir[freq])==1: 
                try: eneend[freq]=float(Edir_min[freq])
                except ValueError:
                    with open(Edir_min[freq]+'_OSZICAR','r') as emfile:
                        line = emfile.readlines()
                    pt = -1
                    while not 'E0' in line[pt]: pt = pt - 1
                    eneend[freq] = float(self.getinfo(line[pt])[4])
        for freq in Edir.keys():
            if len(Edir[freq])==1: enebarr[freq] = Edir[freq][0]
            else: enebarr[freq] = enesaddle[freq] - eneend[freq] 
        return enebarr
    
    def get_v(self,vdir,vdir_num,vdir_denom):
        """Obtaining the attempt frequencies.
            Args:
                vdir <list of dict>: vdir={'v0':*,'v1':*,...}, * is the list like [<local_min>,<saddle>]
                vdir_num <dict>: vdir_num={'E0':<local_min>,...}
                Edir_denum <dict>: Edir_denom={'E0':<saddle>,...}
        """    
        v = dict()
        v_num = dict()
        v_denom = dict()
        for freq in vdir.keys():
            if len(vdir[freq])==1: v[freq] = vdir[freq][0]
            else:
                if os.path.isfile(vdir_num[freq]+'_FREQ') and os.path.isfile(vdir_num[freq]+'_FREQ'): # Reading data from FREQ files
                    with open(vdir_num[freq]+'_FREQ','r') as fn:
                        fnline=fn.readline()
                    nthzlist=self.getinfo(fnline)
                    i=0
                    num_num=0
                    v_num[freq]=1.0
                    while i<len(nthzlist):
                        if abs(float(nthzlist[i]))>0.1: 
                            v_num[freq]=v_num[freq]*float(nthzlist[i])
                            num_num+=1
                        i+=1
                    with open(vdir_denom[freq]+'_FREQ','r') as fd:
                        fdline=fd.readline()
                    dthzlist=self.getinfo(fdline)
                    j=0
                    denom_num=0
                    v_denom[freq]=1.0
                    while j<len(dthzlist):           
                        if abs(float(dthzlist[j]))>0.1: 
                            v_denom[freq]=v_denom[freq]*float(dthzlist[j])
                            denom_num+=1
                        j+=1
                    if not num_num==denom_num+1: 
                        print 'WARNING: Numbers of non-zero frequencies at local minimum and at saddle point do not match! Please double-check the FREQ files of %s and %s'%(vdir_num[freq],vdir_denom[freq])

                else: # Reading data from OUTCAR files
                    nthzlist = fileutil.grepme(vdir_num[freq]+'_OUTCAR', "2PiTHz")
                    num_num=len(nthzlist)
                    dthzlist = fileutil.grepme(vdir_denom[freq]+'_OUTCAR', "2PiTHz")
                    denom_num=len(dthzlist)
                    im_num=0
                    im_denom=0
                    if denom_num==num_num:
                        for i in range(num_num):
                            if 'f/i' in dthzlist[i]: im_denom+=1
                            if 'f/i' in nthzlist[i]: im_num+=1
                        if not im_num==0: 
                            print 'WARNING: Imaginary frequency found in the local minimum! Please check the OUTCAR of %s!'%vdir_num[freq]
                        elif im_denom==0:
                            print 'WARNING: No imaginary frequency found in the saddle point! Please check the OUTCAR of %s!'%vdir_denom[freq]
                            print 'ATTEMPTING TO CORRECT:'
                            numlist=list()
                            for cct in range(num_num):
                                thzfloat = float(dthzlist[cct].split()[3])
                                numlist.append(thzfloat)
                            dthz_std = np.std(numlist)
                            dthz_mean = np.mean(numlist)
                            rangemax = dthz_mean + dthz_std
                            rangemin = dthz_mean - dthz_std
                            removeidx=list()
                            for cct in range(num_num):
                                print numlist[cct]
                                if (numlist[cct] < rangemin):
                                    removeidx.append(cct)
                            if len(removeidx) == 1:
                                print "Removing ", dthzlist[removeidx[0]]
                                dthzlist.pop(removeidx[0])
                                denom_num = denom_num - 1
                                im_denom = im_denom + 1
                            else:
                                print "Could not determine a correction."
                        elif im_denom>1:
                            print 'WARNING: More than 1 imaginary frequencies found in the saddle point! Please check the OUTCAR of %s!'%vdir_denom[freq]
                    else:
                        print 'WARNING: Numbers of frequencies at local minimum and at saddle point are not equal! Please check the OUTCAR of %s and %s!'%(vdir_num[freq],vdir_denom[freq])
                    v_num[freq]=v_denom[freq]=1.0
                    for i in range(num_num):
                        if not 'f/i' in nthzlist[i]:
                            v_num[freq]*=float(self.getinfo(nthzlist[i])[3])
                    for i in range(denom_num):    
                        if not 'f/i' in dthzlist[i]:
                            v_denom[freq]*=float(self.getinfo(dthzlist[i])[3])
                v[freq]=v_num[freq]/v_denom[freq]
        return v


class DiffCoeff(ParsingInputFiles):   
    """Calculating diffusion coefficient.
        Attributes:
            self.init_values <dict of dict>: Initializing values needed to calculate D
            self.calculatingD <void>: printing system information and calculating D in desired frequency model.
        D~1/T data file - Diffusivity.txt is generated.
    """     

    def init_values(self):
        """Evaluating the values that are required in the D calculation, including 
        frequency type, lattice constant, energy barriers, attempt frequencies and HVf (and HB).
        """
        values = dict()      
        vdir_num = dict()
        vdir_denom = dict()
        Edir_saddle = dict()
        Edir_min = dict()
        vdir = self.get_item_name('v')        
        Edir = self.get_item_name('E')
        Hdir = self.get_item_name('H')
        latt = self.get_item_name('lattice')
        model = self.get_item_name('type')
        plotdisplay = self.get_item_name('plotdisplay')
        if len(plotdisplay) == 0: #empty dictionary
            plotdisplay=""
        values['plotdisplay'] = plotdisplay
        numatom = self.get_latt(model,latt)['No.']
        values['a'] = self.get_latt(model,latt)['a']
        values['type'] = model
        if model==8: values['c'] = self.get_latt(model,latt)['c']

        for i in vdir.keys(): 
            if not len(vdir[i])==1:
                try: vdir_num[i]=float(vdir[i][0])
                except ValueError: vdir_num[i] = vdir[i][0]
                try: vdir_denom[i]=float(vdir[i][1])
                except ValueError: vdir_denom[i] = vdir[i][1]
        for i in Edir.keys():
            if not len(Edir[i])==1:
                try: Edir_saddle[i]=float(Edir[i][1])
                except ValueError: Edir_saddle[i] = Edir[i][1]
                try: Edir_min[i] = float(Edir[i][0])
                except ValueError: Edir_min[i] = Edir[i][0]
        values['enebarr'] = self.get_barrier(Edir,Edir_saddle,Edir_min)
        values['v'] = self.get_v(vdir,vdir_num,vdir_denom)
        values['HVf'] = self.get_HB_and_HVf(Hdir,numatom,'HVf')
        if model==8:
            values['HB'] = self.get_HB_and_HVf(Hdir,numatom,'HB')
        return values
    
    def calculatingD(self):
        """Printing system information and calculating D accordingly.
        D~1/T relation is printed and stored in the data file of Diffusivity.txt.
        """
        kB = 11604.52
        jfreq = dict()
        values = self.init_values()
        model = values['type']
        enebarr = values['enebarr']
        v = values['v']
        a = values['a'] 
        plotdisplay = values['plotdisplay']
        HVf = values['HVf'] 
        if model==8: 
            c = values['c'] 
            HB = values['HB'] 
        # print all system information
        if model==5:
            print "FCC Five-Frequency Dilute Diffusion Model"
            print "FCC lattice constant [Angstrom]: {0:.4f}".format(a*10**8)
            print "Energy Barriers [eV]:       E0: {E0:.4f}  E1: {E1:.4f}  E2: {E2:.4f}  E3: {E3:.4f}  E4: {E4:.4f}".format(**enebarr)
            print "Attempt Frequencies [THz]:  v0: {v0:.4f}  v1: {v1:.4f}  v2: {v2:.4f}  v3: {v3:.4f}  v4: {v4:.4f}".format(**v)
            print "Vacancy Formation Energy [eV]: {0:.4f}\n".format(HVf)
            print ""
        if model==8:
            print "HCP Eight-Frequency Dilute Diffusion Model"
            print "HCP basal lattice constant  [Angstrom]: {0:.4f}".format(a*10**8)
            print "HCP c-axis lattice constant [Angstrom]: {0:.4f}".format(c*10**8)
            print "Energy Barriers [eV]:       Ea: {Ea:.4f}  Eb: {Eb:.4f}  Ec: {Ec:.4f}  EX: {EX:.4f}  E'a: {Eap:.4f}  E'b: {Ebp:.4f}  E'c: {Ecp:.4f}  E'X: {EXp:.4f}".format(**enebarr)
            print "Attempt Frequencies [THz]:  va: {va:.4f}  vb: {vb:.4f}  vc: {vc:.4f}  vX: {vX:.4f}  v'a: {vap:.4f}  v'b: {vbp:.4f}  v'c: {vcp:.4f}  v'X: {vXp:.4f}".format(**v)
            print "Vacancy Formation Energy [eV]: {0:.4f}".format(HVf)
            print "Vacancy-Solute Binding Energy [eV]:  {0:.4f}\n".format(HB)
        if model==9:
            print "BCC Nine-Frequency Dilute Diffusion Model"
            print "BCC lattice constant [Angstrom]: {0:.4f}".format(a*10**8)
            print "Energy Barriers [eV]:       E0: {E0:.4f}      E2: {E2:.4f}      E3: {E3:.4f}  E4: {E4:.4f}  E'3: {E3p:.4f}  E'4: {E4p:.4f}".format(**enebarr)
            print "Energy Barriers [eV]:       E''3: {E3pp:.4f}  E''4: {E4pp:.4f}  E5: {E5:.4f}  E6: {E6:.4f}".format(**enebarr)
            print "Attempt Frequencies [THz]:  v0: {v0:.4f}      v2: {v2:.4f}      v3: {v3:.4f}  v4: {v4:.4f}  v'3: {v3p:.4f}  v'4: {v4p:.4f}".format(**v)
            print "Attempt Frequencies [THz]:  v''3: {v3pp:.4f}  v''4: {v4pp:.4f}  v5: {v5:.4f}  v6: {v6:.4f}".format(**v)
            print "Vacancy Formation Energy [eV]: {0:.4f}\n".format(HVf)
            print ""
        if model==14:
            print "FCC 14-Frequency Concentrated Diffusion Model"
            print "FCC lattice constant [Angstrom]: {0:.4f}".format(a*10**8)
            print "Energy Barriers [eV]:       E0: {E0:.4f}    E1: {E1:.4f}    E2: {E2:.4f}    E3: {E3:.4f}    E4: {E4:.4f}".format(**enebarr)
            print "Energy Barriers [eV]:       E11: {E11:.4f}  E12: {E12:.4f}  E13: {E13:.4f}  E14: {E14:.4f}".format(**enebarr)
            print "Energy Barriers [eV]:       E23: {E23:.4f}  E24: {E24:.4f}  E33: {E33:.4f}  E34: {E34:.4f}  E44: {E44:.4f}".format(**enebarr)
            print "Attempt Frequencies [THz]:  v0: {v0:.4f}    v1: {v1:.4f}    v2: {v2:.4f}    v3: {v3:.4f}    v4: {v4:.4f}".format(**v)
            print "Attempt Frequencies [THz]:  v11: {v11:.4f}  v12: {v12:.4f}  v13: {v13:.4f}  v14: {v14:.4f}".format(**v)
            print "Attempt Frequencies [THz]:  v23: {v23:.4f}  v24: {v24:.4f}  v33: {v33:.4f}  v34: {v34:.4f}  v44: {v44:.4f}".format(**v)
            print "Vacancy Formation Energy [eV]: {0:.4f}\n".format(HVf)
            print ""
        try:
            temp = self.get_item_name('temp')['temp']
            tempstart = temp[0]
            tempstep = temp[1]
            tempend = temp[2]
        except: tempstart = 0.0; tempstep = 0.1; tempend = 2.0 # default temperature range
        if model==5:
            with open('Diffusivity.txt','w+') as fp:
                fp.write('1000/T[K^(-1)]    D_solute[cm^2/s]\n')
            print '1000/T[K^(-1)]    D_solute[cm^2/s]'
            D = []
            f0 = 0.7815
            t = tempstart
            i = 0
            while t<tempend+tempstep:
                for freq in ['v0','v1','v2','v3','v4']:
                    if t==0.0: jfreq[freq] = v[freq]*10**12
                    else: 
                        T=1000./t
                        jfreq[freq] = v[freq]* 10**12 * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
                alpha = jfreq['v4']/jfreq['v0']
                F_num = 10*np.power(alpha,4) + 180.5*np.power(alpha,3) + 927*np.power(alpha,2) + 1341*alpha
                F_denom = 2*np.power(alpha,4) + 40.2*np.power(alpha,3) + 254*np.power(alpha,2) + 597*alpha + 435
                FX = 1-(1.0/7.0)*(F_num/F_denom)
                f2_num = 1+3.5*FX *(jfreq['v3']/jfreq['v1'])
                f2_denom = 1+(jfreq['v2']/jfreq['v1']) + 3.5*FX*(jfreq['v3']/jfreq['v1'])
                f2 = f2_num / f2_denom       
                if t==0.0: Vacconc=1.0
                else: Vacconc = np.exp(-kB*HVf/T)       
                Dself = f0*Vacconc*a**2*jfreq['v0']
                D.append(Dself*f2*jfreq['v2']*jfreq['v4']/f0/jfreq['v0']/jfreq['v3'])
                with open('Diffusivity.txt','a+') as fp:
                    fp.write('%f   %E\n'%(t,D[i]))
                print '%f   %E'%(t,D[i])
                t = t + tempstep    
                i = i + 1
            if plotdisplay.lower() == "none": #No plot display
                pass
            else:
                import matplotlib
                if len(plotdisplay) == 0:
                    pass
                else:
                    matplotlib.use(plotdisplay)
                import matplotlib.pyplot as plt
                tvals = np.linspace(tempstart, tempend, i)
                plt.plot(tvals, np.log10(D), '.')
                plt.xlabel('$10^3/T$ (K$^{-1}$)')
                plt.ylabel('log$_{10}$$D$ (cm$^2$/s)')
                plt.savefig('Diffusivity.png')
                plt.show()
        if model==9:
            with open('Diffusivity.txt','w+') as fp:
                fp.write('1000/T[K^(-1)]    D_solute[cm^2/s]\n')
            print '1000/T[K^(-1)]    D_solute[cm^2/s]'
            D = []
            f0 = 0.7272
            FF = 0.512
            t = tempstart
            i = 0
            while t<tempend+tempstep:
                for freq in ['v0','v2','v3','v4','v3p','v4p','v3pp','v4pp','v5','v6']:
                    if t==0.0: jfreq[freq] = v[freq]*10**12
                    else: 
                        T=1000./t
                        jfreq[freq] = v[freq]* 10**12 * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
                tt = -jfreq['v2']/( jfreq['v2'] + 3*jfreq['v3'] + 3*jfreq['v3p'] + jfreq['v3pp'] - (jfreq['v3']*jfreq['v4'])/(jfreq['v4']+FF*jfreq['v5']) - (2*jfreq['v3p']*jfreq['v4p'])/(jfreq['v4p']+3*FF*jfreq['v0']) - (jfreq['v3pp']*jfreq['v4pp'])/(jfreq['v4pp']+7*FF*jfreq['v0']) );
                f2 = (1 + tt) / (1 - tt)       
                if t==0.0: Vacconc=1.0
                else: Vacconc = np.exp(-kB*HVf/T)       
                Dself = f0*Vacconc*a**2*jfreq['v0']
                D.append(a**2*Vacconc*f2*jfreq['v2']*jfreq['v4p']/jfreq['v3p'])
                with open('Diffusivity.txt','a+') as fp:
                    fp.write('%f   %E\n'%(t,D[i]))
                print '%f   %E'%(t,D[i])
                t = t + tempstep    
                i = i + 1
            if plotdisplay.lower() == "none": #No plot display
                pass
            else:
                import matplotlib
                if len(plotdisplay) == 0:
                    pass
                else:
                    matplotlib.use(plotdisplay)
                import matplotlib.pyplot as plt
                tvals = np.linspace(tempstart, tempend, i)
                plt.plot(tvals, np.log10(D), '.')
                plt.xlabel('$10^3/T$ (K$^{-1}$)')
                plt.ylabel('log$_{10}$$D$ (cm$^2$/s)')
                plt.savefig('Diffusivity.png')
                plt.show()
        if model==14:
            with open('Diffusivity.txt','w+') as fp:
                fp.write('D_host(c) = D_host(0) * [ 1 + b1*c + b2*c^2 + ...]\n')
                fp.write('D_solute(c) = D_solute(0) * [ 1 + B1*c + ...]\n')
                fp.write('c = solute concentration [atomic fraction]\n\n')
                fp.write('1000/T[K^(-1)]    D_solute[cm^2/s]    b1        b2        B1\n')
            print 'D_host(c) = D_host(0) * [ 1 + b1*c + b2*c^2 + ...]'
            print 'D_solute(c) = D_solute(0) * [ 1 + B1*c + ...]'
            print 'c = solute concentration [atomic fraction]\n'
            print '1000/T[K^(-1)]    D_solute[cm^2/s]    b1        b2        B1'
            D = []
            b1 = []
            b2 = []
            B1 = []
            f0 = 0.7815
            t = tempstart
            i = 0
            while t<tempend+tempstep:
                for freq in ['v0','v1','v2','v3','v4','v11','v12','v13','v14','v23','v24','v33','v34','v44']:
                    if t==0.0: jfreq[freq] = v[freq]*10**12
                    else: 
                        T=1000./t
                        jfreq[freq] = v[freq]* 10**12 * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
                alpha = jfreq['v4']/jfreq['v0']
                F_num = 10*np.power(alpha,4) + 180.5*np.power(alpha,3) + 927*np.power(alpha,2) + 1341*alpha
                F_denom = 2*np.power(alpha,4) + 40.2*np.power(alpha,3) + 254*np.power(alpha,2) + 597*alpha + 435
                FX = 1-(1.0/7.0)*(F_num/F_denom)
                f2_num = 1+3.5*FX *(jfreq['v3']/jfreq['v1'])
                f2_denom = 1+(jfreq['v2']/jfreq['v1']) + 3.5*FX*(jfreq['v3']/jfreq['v1'])
                f2 = f2_num / f2_denom       
                if t==0.0: Vacconc=1.0
                else: Vacconc = np.exp(-kB*HVf/T)
                Dself = f0*Vacconc*a**2*jfreq['v0']
                D.append(Dself*f2*jfreq['v2']*jfreq['v4']/f0/jfreq['v0']/jfreq['v3'])
                b1.append(-18 + (jfreq['v4']/jfreq['v3'])*((4*jfreq['v1'] + 14*jfreq['v3'])/jfreq['v0']))
                b2.append((jfreq['v23']/jfreq['v24'])*(jfreq['v4']/jfreq['v3'])*((jfreq['v4']/jfreq['v3'])*(jfreq['v4']/jfreq['v3'])*(8*jfreq['v13']+2*jfreq['v11']+10*jfreq['v33']-16*jfreq['v1']-28*jfreq['v3'])/jfreq['v0'] + (jfreq['v4']/jfreq['v3'])*(8*(jfreq['v14']-jfreq['v1'])+2*(jfreq['v34']-jfreq['v3']))/jfreq['v0'] +(10*jfreq['v44']-37*jfreq['v4'])/jfreq['v0'] + 63) + (jfreq['v4']/jfreq['v3'])*(jfreq['v4']/jfreq['v3'])*(4*jfreq['v11']+20*jfreq['v13']+11*jfreq['v33']-28*jfreq['v1']-49*jfreq['v3'])/jfreq['v0'] + (jfreq['v4']/jfreq['v3'])*(20*(jfreq['v14']-jfreq['v1'])+47*(jfreq['v34']-jfreq['v3']))/jfreq['v0'] + 11*jfreq['v44']/jfreq['v0'] - 173*jfreq['v4']/jfreq['v0'] + 192)
                B1.append(-6 - 12*((jfreq['v23']/jfreq['v24'])*(jfreq['v4']/jfreq['v3'])) + (jfreq['v4']/jfreq['v3'])*(4*(jfreq['v21']/jfreq['v2'])*(jfreq['v23']/jfreq['v24'])*(jfreq['v4']/jfreq['v3']) + 14*(jfreq['v23']/jfreq['v2'])))
                with open('Diffusivity.txt','a+') as fp:
                    fp.write('%f   %E    %f    %f    %f\n'%(t,D[i],b1[i],b2[i],B1[i]))
                print '%f   %E    %f    %f    %f'%(t,D[i],b1[i],b2[i],B1[i])
                t = t + tempstep    
                i = i + 1
            if plotdisplay.lower() == "none": #No plot display
                pass
            else:
                import matplotlib
                if len(plotdisplay) == 0:
                    pass
                else:
                    matplotlib.use(plotdisplay)
                import matplotlib.pyplot as plt
                tvals = np.linspace(tempstart, tempend, i)
                plt.plot(tvals, np.log10(D), '.')
                plt.xlabel('$10^3/T$ (K$^{-1}$)')
                plt.ylabel('log$_{10}$$D$ (cm$^2$/s)')
                plt.savefig('Diffusivity.png')
                plt.show()
        elif model==8:
            with open('Diffusivity.txt','w+') as fp:
                fp.write('1000/T[K^(-1)]    D_basal[cm^2/s]    D_c-axis[cm^2/s]\n')
            print '1000/T[K^(-1)]    D_basal[cm^2/s]    D_c-axis[cm^2/s]'
            Dperp = []
            Dpara = []
            F = 0.736
            lambda1 = a
            lambda2b = a/(3**0.5)
            t = tempstart
            i = 0
            while t<tempend+tempstep:
                for freq in ['va','vb','vc','vX','vap','vbp','vcp','vXp']:
                    if t==0.0: jfreq[freq] = v[freq]*10**12
                    else:
                        T=1000./t
                        jfreq[freq] = v[freq] * 10**12 * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
                f1z = (2*jfreq['vap']+7*F*jfreq['vcp'])/(2*jfreq['vap']+7*F*jfreq['vcp']+2*jfreq['vXp'])
                S1b_num = (3**0.5)*jfreq['vX']/jfreq['va']*lambda1+jfreq['vXp']/jfreq['vap']*lambda2b*(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
                S1b_denom = 3 - (2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])*(2+3*jfreq['vbp']/jfreq['vap']+7*F*jfreq['vcp']/jfreq['vap']+2*jfreq['vXp']/jfreq['vap']-1./(2+3*jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+jfreq['vX']/jfreq['va']))
                S1b = S1b_num / S1b_denom
                S2x=((3**0.5)*S1b-jfreq['vX']/jfreq['va']*lambda1)/(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
                f1b=1+2*S1b/lambda2b
                f2x=1+2*S2x/lambda1
                if t==0.0: C2=1.0
                else: C2=np.exp(-kB*(HVf+HB)/T)
                Dperp.append(a**2*(3*f2x*jfreq['vX']+f1b*jfreq['vXp'])*C2/2)
                Dpara.append(c**2*f1z*jfreq['vXp']*C2*0.75)
                with open('Diffusivity.txt','a+') as fp:
                    fp.write('%f   %E   %E\n'%(t,Dperp[i],Dpara[i]))
                print '%f   %E   %E'%(t,Dperp[i],Dpara[i])
                t = t + tempstep
                i = i + 1
            if plotdisplay.lower() == "none":
                pass
            else:
                import matplotlib
                if len(plotdisplay) == 0:
                    pass
                else:
                    matplotlib.use(plotdisplay)
                import matplotlib.pyplot as plt
                tvals = np.linspace(tempstart, tempend, i)
                plt.plot(tvals, np.log10(Dperp), '.', tvals, np.log10(Dpara), '+')
                plt.xlabel('$10^3/T$ (K$^{-1}$)')
                plt.ylabel('log$_{10}$$D$ (cm$^2$/s)')
                plt.savefig('Diffusivity.png')
                plt.show()

