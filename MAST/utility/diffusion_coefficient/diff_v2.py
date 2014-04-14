#!/usr/bin/env python
import sys, getopt, os
import numpy as np
import pymatgen as mg
from MAST.utility import fileutil
def getinfo(line):
    line=line.strip('\n')
    data=line.split(' ')
    while 1:
        try: data.remove('')
        except ValueError: break
    return data

jfreq = {}
v = {}
v_num = {}
v_denom = {}
enebarr = {}
enesaddle = {}
eneend = {}
vdir_num = {}
vdir_denom = {}
Edir_neb = {}
Edir_def = {}
def get_item_name(inp,keyword):
    item_name = {}
    content = open(inp,'r').readlines()
    Elist=['E0','E1','E2','E3','E4','Ea','Eb','Ec','EX','Eap','Ebp','Ecp','EXp']
    vlist=['v0','v1','v2','v3','v4','va','vb','vc','vX','vap','vbp','vcp','vXp']
    Hlist=['HB','HVf']
    for i in range(len(content)):
        line = getinfo(content[i])
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
                    item_name['type'] = 5
                elif line[0]=='type' and line[1]=='hcp_8freq':  
                    item_name['type'] = 8
            elif keyword=='lattice' and 'lattice' in line[0]:
                item_name['lattice'] = line[1]
    return item_name

def get_latt(lattice):
    data={'a':0,'c':0,'No.':0}
    os.system('cp '+lattice+'_POSCAR POSCAR')
    struct = mg.read_structure('POSCAR')
    os.system('rm POSCAR')
    reduced = mg.symmetry.finder.SymmetryFinder(struct,0.001).get_primitive_standard_structure()   
    data['No.'] = len(struct)
    if types['type']==5:
        data['a'] = reduced.lattice.abc[0]*np.sqrt(2)*10**(-8)
    if types['type']==8:
        data['a'] = reduced.lattice.abc[0]*10**(-8)
        data['c'] = reduced.lattice.abc[2]*10**(-8)
    return data
    
argv = sys.argv[1:]
opts, args = getopt.getopt(argv,"i")
inp = args[0]
vdir = get_item_name(inp,'v')
types = get_item_name(inp,'type')
Edir = get_item_name(inp,'E')
Hdir = get_item_name(inp,'H')
lattice = get_item_name(inp,'lattice')['lattice']
numatom = get_latt(lattice)['No.']
# need to get MAST to pass on POSCAR as well as OSZICAR for the following to be able to change numatom behavior properly
#perfect = get_item_name(inp,'H')['HVf']['perfect']
#vacancy = get_item_name(inp,'H')['HVf']['vac']
#numatom_perfect = get_latt(perfect)['No.']
#numatom_vacancy = get_latt(vacancy)['No.']
a = get_latt(lattice)['a']
if types['type']==8: c = get_latt(lattice)['c']

for i in vdir.keys(): 
    if not len(vdir[i])==1:
        try: vdir_num[i]=float(vdir[i][0])
        except ValueError: vdir_num[i] = vdir[i][0]
        try: vdir_denom[i]=float(vdir[i][1])
        except ValueError: vdir_denom[i] = vdir[i][1]
for i in Edir.keys():
    if not len(Edir[i])==1:
        try: Edir_neb[i]=float(Edir[i][1])
        except ValueError: Edir_neb[i] = Edir[i][1]
        try: Edir_def[i] = float(Edir[i][0])
        except ValueError: Edir_def[i] = Edir[i][0]

def get_HB_and_HVf(Hdir,keyword):
    try: 
        float(Hdir[keyword])
        return Hdir[keyword]
    except TypeError:
        dir = {}
        ene = {}
        for key in Hdir[keyword].keys():
            try: 
                ene[key]=float(Hdir[keyword][key])
            except ValueError: 
                dir[key] = Hdir[keyword][key]
                line = open(dir[key]+'_OSZICAR','r').readlines()      
                pt = -1   
                while not 'E0' in line[pt]: pt = pt - 1
                ene[key] = float(getinfo(line[pt])[4])
        if keyword=='HVf': #HHW numatom behavior need to be changed for HVf calculation
            return ene['vac'] - (numatom - 1)*ene['perfect']/numatom
            # return ene['vac'] - (numatom_vacancy)*ene['perfect']/numatom_perfect
        elif keyword=='HB': #HHW additional checks for numatom will be needed for HB calculation (though not as likely)
            return ene['perfect'] + ene['vac-sub'] - ene['sub'] - ene['vac']

def get_saddle(Edir,Edir_neb):
    for freq in Edir.keys():
        if not len(Edir[freq])==1: 
            try: enesaddle[freq]=float(Edir_neb[freq])
            except ValueError:
                line = open(Edir_neb[freq]+'_OSZICAR','r').readlines()
                pt = -1
                while not 'E0' in line[pt]: pt = pt - 1
                enesaddle[freq] = float(getinfo(line[pt])[4])
    return enesaddle

def get_end(Edir,Edir_def):
    for freq in Edir.keys():
        if not len(Edir[freq])==1: 
            try: eneend[freq]=float(Edir_def[freq])
            except ValueError:
                line = open(Edir_def[freq]+'_OSZICAR','r').readlines()
                pt = -1
                while not 'E0' in line[pt]: pt = pt - 1
                eneend[freq] = float(getinfo(line[pt])[4])
    return eneend

#HHW need to add check for number of frequencies found
def get_v(vdir,vdir_num,vdir_denom):
    for freq in vdir.keys():
        if len(vdir[freq])==1: v[freq] = vdir[freq][0]*10**12
        else:
            try: v_num[freq]=float(vdir_num[freq])
            except ValueError:
                if os.path.isfile(vdir_num[freq]+'_FREQ'):
                    fn=open(vdir_num[freq]+'_FREQ','r')
                    freqn=getinfo(fn.readline())
                    i=0
                    v_num[freq]=1.0
                    while i<len(freqn):
                        if not float(freqn[i])==0.0: v_num[freq]=v_num[freq]*float(freqn[i])
                        i+=1
                    fn.close()
                else: #TTM add OUTCAR direct reading block
                    fname = vdir_num[freq]+'_OUTCAR'
                    thzlist = fileutil.grepme(fname, "THz")
                    v_num[freq]=1.0
                    for idx in range(0, len(thzlist)):
                        thzline = thzlist[idx]
                        thzsplit = thzline.split()
                        label = thzsplit[1]
                        if not 'i' in label:
                            v_num[freq] = v_num[freq]*float(thzsplit[3])
            try: v_denom[freq]=float(vdir_denom[freq])
            except ValueError:
                if os.path.isfile(vdir_num[freq]+'_FREQ'):
                    fd=open(vdir_denom[freq]+'_FREQ','r')
                    freqd=getinfo(fd.readline())
                    i=0
                    v_denom[freq]=1.0
                    while i<len(freqd):           
                        if not float(freqd[i])==0.0: v_denom[freq]=v_denom[freq]*float(freqd[i])
                        i+=1
                    fd.close()
                else: #TTM add OUTCAR direct reading block
                    fname = vdir_denom[freq]+'_OUTCAR'
                    thzlist = fileutil.grepme(fname, "THz")
                    v_denom[freq]=1.0
                    for idx in range(0, len(thzlist)):
                        thzline = thzlist[idx]
                        thzsplit = thzline.split()
                        label = thzsplit[1]
                        if not 'f/i' in label:
                            v_denom[freq] = v_denom[freq]*float(thzsplit[3])
            v[freq]=v_num[freq]/v_denom[freq]*10**12
    return v


kB = 11604.52
for freq in Edir.keys():
    if len(Edir[freq])==1: enebarr[freq] = Edir[freq][0]
    else: enebarr[freq] = get_saddle(Edir,Edir_neb)[freq] - get_end(Edir,Edir_def)[freq]
v = get_v(vdir,vdir_num,vdir_denom)
v_THz = v.copy()  # attempt freq in THz for screen output
for freq in vdir.keys():
    v_THz[freq] = v_THz[freq]*10**(-12)
HVf = get_HB_and_HVf(Hdir,'HVf')
if types['type']==8:
    HB = get_HB_and_HVf(Hdir,'HB')


# print all system information
if types['type']==5:
    print "FCC Five-Frequency Diffusion Model"
    print "FCC lattice constant [angstrom]: {0:.4f}".format(a*10**8)
    print "Energy Barriers [eV]:       E0: {E0:.4f}  E1: {E1:.4f}  E2: {E2:.4f}  E3: {E3:.4f}  E4: {E4:.4f}".format(**enebarr)
    print "Attempt Frequencies [THz]:  v0: {v0:.4f}  v1: {v1:.4f}  v2: {v2:.4f}  v3: {v3:.4f}  v4: {v4:.4f}".format(**v_THz)
    print "Vacancy Formation Energy [eV]: {0:.4f}\n".format(HVf)
    print ""
if types['type']==8:
    print "HCP Eight-Frequency Diffusion Model"
    print "HCP basal lattice constant  [angstrom]: {0:.4f}".format(a*10**8)
    print "HCP c-axis lattice constant [angstrom]: {0:.4f}".format(c*10**8)
    print "Energy Barriers [eV]:       Ea: {Ea:.4f}  Eb: {Eb:.4f}  Ec: {Ec:.4f}  EX: {EX:.4f}  E'a: {Eap:.4f}  E'b: {Ebp:.4f}  E'c: {Ecp:.4f}  E'X: {EXp:.4f}".format(**enebarr)
    print "Attempt Frequencies [THz]:  va: {va:.4f}  vb: {vb:.4f}  vc: {vc:.4f}  vX: {vX:.4f}  v'a: {vap:.4f}  v'b: {vbp:.4f}  v'c: {vcp:.4f}  v'X: {vXp:.4f}".format(**v_THz)
    print "Vacancy Formation Energy [eV]: {0:.4f}".format(HVf)
    print "Vacancy-Solute Binding Energy [eV]:  {0:.4f}\n".format(HB)


try:
    temp = get_item_name(inp,'temp')['temp']
    t0 = temp[0]
    dt = temp[1]
    t1 = temp[2]
except: t0 = 0.0; dt = 0.1; t1 = 2.0
if types['type']==5:
    fp = open('Diffusivity.txt','w+')
    fp.write('1000/T(K^(-1))    D(cm^2/s)\n')
    print '1000/T(K^(-1))    D(cm^2/s)'
    D = []
    f0 = 0.7815
    t = t0
    i = 0
    while t<t1+dt:
        for freq in vdir.keys():
            if t==0.0: jfreq[freq] = v[freq]
            else: 
                T=1000./t
                jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
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
        fp.write('%f  %E\n'%(t,D[i]))
        print '%f  %E'%(t,D[i])
        t = t + dt    
        i = i + 1
    
elif types['type']==8:
    fp = open('Diffusivity.txt','w+')
    fp.write('1000/T(K^(-1))    D_basal(cm^2/s)    D_c-axis(cm^2/s)\n')
    print '1000/T(K^(-1))    D_basal(cm^2/s)    D_c-axis(cm^2/s)'
    Dperp = []
    Dpara = []
    F = 0.736
    lambda1 = a
    lambda2b = a/(3**0.5)
    t = t0
    i = 0
    while t<t1+dt:
        for freq in vdir.keys():
            if t==0.0: jfreq[freq] = v[freq]
            else:
                T=1000./t
                jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
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
        fp.write('%f  %E  %E\n'%(t,Dperp[i],Dpara[i]))
        print '%f  %E  %E'%(t,Dperp[i],Dpara[i])
        t = t + dt
        i = i + 1
    


