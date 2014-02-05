#!/usr/bin/env python
import sys, getopt, os
import numpy as np
import matplotlib.pyplot as plt
import pymatgen as mg

def getinfo(line):
    line=line.strip('\n')
    data=line.split(' ')
    while 1:
        try: data.remove('')
        except: break
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
def get_freq_name(inp,keyword):
    freq_name = {}
    content = []
    flag = 0
    fp = open(inp,'r')
    for line in fp:
        #TTM 20140130 separate input file; no $freq/$end
        #if '$freq' in line:
        #    flag = 1
        #if '$end' in line:
        #    flag = 0
        #if flag==1:
        #    content.append(line)
        content.append(line)
    #content.pop(0) #TTM 20140130 separate input file, no $freq
    fp.close()
    #TTM 20140130 mixed keywords; explicitly allow some
    Elist=['E0','E1','E2','E3','E4','Ea','Eb','Ec','EX','Eap','Ebp','Ecp','EXp']
    vlist=['v0','v1','v2','v3','v4','va','vb','vc','vX','vap','vbp','vcp','vXp']
    Hlist=['HB','HVf']
    for i in range(len(content)):
        line = getinfo(content[i])
        if not len(line)==0:
            if (keyword=='E' and line[0] in Elist) or (keyword=='v' and line[0] in vlist):
                freq_name[line[0]] = []
                if len(line)==2:
                    try:
                        freq_name[line[0]].append(float(line[1]))
                    except ValueError:
                        pass
                else:
                    for index in range(1,len(line)):
                        freq_name[line[0]].append(line[index])
            elif (keyword=='H' and line[0] in Hlist):
                if line[0] == 'HB':
                    if len(line) > 2:
                        freq_name['HB'] = dict()
                        freq_name['HB']['perfect'] = line[1]
                        freq_name['HB']['sub'] = line[2]
                        freq_name['HB']['vac-sub'] = line[3]
                        freq_name['HB']['vac'] = line[4]
                    else:
                        freq_name['HB'] = float(line[1])
                elif line[0] == 'HVf':
                    if len(line) > 2:
                        freq_name['HVf'] = dict()
                        freq_name['HVf']['perfect'] = line[1]
                        freq_name['HVf']['vac'] = line[2]
                    else:
                        freq_name['HB'] = float(line[1])
            elif keyword=='temp' and 'temp' in line[0]:
                freq_name['temp'] = [float(line[1]),float(line[2]),float(line[3])]
            elif keyword=='type' and 'type' in line[0]:
                if line[0]=='type' and line[1]=='5' or line[1]=='fcc':
                    freq_name['type'] = 5
                elif line[0]=='type' and line[1]=='8' or line[1]=='hcp':  
                    freq_name['type'] = 8
            elif keyword=='lattice' and 'lattice' in line[0]:
                freq_name['lattice'] = line[1]
    return freq_name
               
def get_freq_sub(dir):
    sub = []
    for folder in os.listdir(dir):
        if 'FREQ' in folder:
            sub.append(folder)
    return 'FREQ'+str(len(sub))

def get_latt():
    data={'a':0,'c':0,'No.':0}
    data['No.'] = len(mg.read_structure(lattice +'/POSCAR'))
    #struct = mg.read_structure(lattice +'/POSCAR_primitive')        #TTM use pymatgen get_primitive_structure
    struct = mg.read_structure(lattice + '/POSCAR').get_primitive_structure()
      
    if round(struct.lattice.angles[0])==round(struct.lattice.angles[1])==round(struct.lattice.angles[2])==90.0:
        data['a'] = struct.lattice.abc[0]*10**(-8)
    for i in range(3):
        if not round(struct.lattice.angles[i])==90:
            data['a'] = struct.lattice.abc[(i+1)%3]*10**(-8)
            data['c'] = struct.lattice.abc[i]*10**(-8)
    return data
argv = sys.argv[1:]
opts, args = getopt.getopt(argv,"i")
inp = os.path.join(os.getcwd(),args[0])
vdir = get_freq_name(inp,'v')
types = get_freq_name(inp,'type')
Edir = get_freq_name(inp,'E')
Hdir = get_freq_name(inp,'H')
lattice = get_freq_name(inp,'lattice')['lattice'] #TTM added 20140129
for i in vdir.keys(): 
    if len(vdir[i])==1: continue
    vdir_denom[i] = os.path.join(os.getcwd(),'phonon_'+vdir[i][1]+'_parse')
    vdir_num[i] = os.path.join(os.getcwd(),'phonon_'+vdir[i][0]+'_parse')
for i in Edir.keys():
    if len(Edir[i])==1: continue
    Edir_neb[i] = os.path.join(os.getcwd(),'neb_'+Edir[i][1]+'_stat')
    Edir_def[i] = os.path.join(os.getcwd(),'defect_'+Edir[i][0]+'_stat')

def get_HB_and_HVf(Hdir,keyword):
    if len(Hdir[keyword])==1: return Hdir[keyword][0]
    else:
        dir = {}
        ene = {}
        for key in Hdir[keyword].keys():
            dir[key] = os.path.join(os.getcwd(),Hdir[keyword][key]+'_stat')         
            for line in open(os.path.join(dir[key],'OSZICAR'),'r'):
                if 'E0' in line:
                    ene[key] = float(getinfo(line)[4])
        if keyword=='HVf':
            return ene['vac'] - (numatom - 1)*ene['perfect']/numatom
        elif keyword=='HB':
            return ene['perfect'] + ene['vac-sub'] - ene['sub'] - ene['vac']

def get_saddle(Edir,Edir_neb):
    for freq in Edir.keys():
        energy = []
        images = 1
        if len(Edir[freq])==1: continue
        for folder in os.listdir(Edir_neb[freq]):
            if '0'+str(images+2) in folder:
                images = images + 1               
        getsaddle=[os.path.join(Edir_neb[freq],'0'+str(j),'OSZICAR') for j in range(images+2)]
        for j in range(images+2):
            free = open(getsaddle[j],'r')
            for line in free:
                if 'E0' in line:
                    energy.append(float(getinfo(line)[4]))
            free.close()
        enesaddle[freq] = max(energy)
    return enesaddle

def get_end(Edir,Edir_def):
    for freq in Edir.keys():
        if len(Edir[freq])==1: continue
        fp = open(os.path.join(Edir_def[freq],'OSZICAR'),'r')
        for line in fp:
            if 'E0' in line:
                eneend[freq] = float(getinfo(line)[4])
        fp.close()
    return eneend

def get_v(vdir,vdir_num,vdir_denom):
    for freq in vdir.keys():
        if len(vdir[freq])==1: v[freq] = vdir[freq][0]*10**12
        else:
            fn=open(os.path.join(vdir_num[freq],get_freq_sub(vdir_num[freq])),'r')
            freqn=getinfo(fn.readline())
            fd=open(os.path.join(vdir_denom[freq],get_freq_sub(vdir_denom[freq])),'r')
            freqd=getinfo(fd.readline())
            n,d=0,0
            v_num[freq],v_denom[freq] = 1.0,1.0
            while n<len(freqn):
                if not float(freqn[n])==0.0: v_num[freq]=v_num[freq]*float(freqn[n])
                n=n+1
            while d<len(freqd):           
                if not float(freqd[d])==0.0: v_denom[freq]=v_denom[freq]*float(freqd[d])
                d=d+1
            v[freq]=v_num[freq]/v_denom[freq]*10**12
    return v


kB = 11604.52
for freq in Edir.keys():
    if len(Edir[freq])==1: enebarr[freq] = Edir[freq][0]
    else: enebarr[freq] = get_saddle(Edir,Edir_neb)[freq] - get_end(Edir,Edir_def)[freq]
v = get_v(vdir,vdir_num,vdir_denom)
numatom = get_latt()['No.']
a = get_latt()['a']
c = get_latt()['c']
HVf = get_HB_and_HVf(Hdir,'HVf')
if types['type']==8:
    HB = get_HB_and_HVf(Hdir,'HB')
fp = open('Diffusivity','w')
try:
    temp = get_freq_name(inp,'temp')['temp']
    t0 = temp[0]
    dt = temp[1]
    t1 = temp[2]
except: t0 = 0.1; dt = 0.1; t1 = 2.0
if types['type']==5:
    fp.write('1000/T(K^(-1))    log D(cm^2/s)\n')
    D = []
    f0 = 0.7815
    t = t0
    i = 0
    while t<t1+dt:
        T=1000./t
        for freq in vdir.keys():
           jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
        alpha = jfreq['v4']/jfreq['v0']
        F_num = 10*np.power(alpha,4) + 180.5*np.power(alpha,3) + 927*np.power(alpha,2) + 1341*alpha
        F_denom = 2*np.power(alpha,4) + 40.2*np.power(alpha,3) + 254*np.power(alpha,2) + 597*alpha + 435
        FX = 1-(1.0/7.0)*(F_num/F_denom)
        f2_num = 1+3.5*FX *(jfreq['v3']/jfreq['v1'])
        f2_denom = 1+(jfreq['v2']/jfreq['v1']) + 3.5*FX*(jfreq['v3']/jfreq['v1'])
        f2 = f2_num / f2_denom       
        Vacconc = np.exp(-kB*HVf/T)       
        Dself = f0*Vacconc*a**2*jfreq['v0']
        D.append(np.log10(Dself*f2*jfreq['v2']*jfreq['v4']/f0/jfreq['v0']/jfreq['v3']))
        fp.write(str(t)+'    '+str(D[i])+'\n')
        t = t + dt    
        i = i + 1
    t = np.linspace(t0,t1,i)    
    plt.plot(t,D,'.')
    plt.xlabel('$10^3/T$(K$^{-1}$)')
    plt.ylabel('log $D$(cm$^2$/s)')
    plt.show()
    
elif types['type']==8:
    fp.write('1000/T(K^(-1))    log D_basal(cm^2/s)    log D_c-axis(cm^2/s)\n')
    Dperp = []
    Dpara = []
    F = 0.736
    lambda1 = a
    lambda2b = a/(3**0.5)
    t = t0
    i = 0
    while t<t1+dt:
        T=1000./t
        for freq in vdir.keys():
           jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*kB/T)
        f1z = (2*jfreq['vap']+7*F*jfreq['vcp'])/(2*jfreq['vap']+7*F*jfreq['vcp']+2*jfreq['vXp'])
        S1b_num = (3**0.5)*jfreq['vX']/jfreq['va']*lambda1+jfreq['vXp']/jfreq['vap']*lambda2b*(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
        S1b_denom = 3 - (2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])*(2+3*jfreq['vbp']/jfreq['vap']+7*F*jfreq['vcp']/jfreq['vap']+2*jfreq['vXp']/jfreq['vap']-1./(2+3*jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+jfreq['vX']/jfreq['va']))
        S1b = S1b_num / S1b_denom
        S2x=((3**0.5)*S1b-jfreq['vX']/jfreq['va']*lambda1)/(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
        f1b=1+2*S1b/lambda2b
        f2x=1+2*S2x/lambda1
        C2=np.exp(-kB*(HVf+HB)/T)
        Dperp.append(np.log10(a**2*(3*f2x*jfreq['vX']+f1b*jfreq['vXp'])*C2/2))
        Dpara.append(np.log10(c**2*f1z*jfreq['vXp']*C2*0.75))
        fp.write(str(t)+'    '+str(Dperp[i])+'    '+str(Dpara[i])+'\n')
        t = t + dt
        i = i + 1
    print t
    t = np.linspace(t0,t1,i)
    plt.plot(t,Dperp,'.',t,Dpara,'+')
    plt.xlabel('$10^3/T$(K$^{-1}$)')
    plt.ylabel('log $D$(cm$^2$/s)')
    plt.show()
