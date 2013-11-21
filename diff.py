import os
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
per_dir = os.path.join(os.getcwd(),'perfect_stat')
def_dir = os.path.join(os.getcwd(),'defect_vac_stat')

def get_freq_name(inp,keyword):
    freq_name = {}
    content = []
    flag = 0
    fp = open(inp,'r')
    for line in fp:
        if '$path' in line:
            flag = 1
        if '$end' in line:
            flag = 0
        if flag==1:
            content.append(line)
    content.pop(0)
    fp.close()
    for i in range(len(content)):
        line = getinfo(content[i])
        if not len(line)==0:
            if keyword=='E' and 'E' in line[0] or keyword=='v' and 'v' in line[0]:
                freq_name[line[0]] = []
                freq_name[line[0]].append(line[1])
                freq_name[line[0]].append(line[2])
    return freq_name
               
def get_freq_sub(dir):
    sub = []
    for folder in os.listdir(dir):
	if 'FREQ' in folder:
	    sub.append(folder)
    return 'FREQ'+str(len(sub))

def get_latt():
    data={'a':0,'c':0,'No.':0}
    data['No.'] = len(mg.read_structure(os.getcwd()+'/POSCAR'))
    struct = mg.read_structure(os.getcwd()+'/POSCAR_primitive')    
      
    if round(struct.lattice.angles[0])==round(struct.lattice.angles[1])==round(struct.lattice.angles[2])==90.0:
        data['a'] = struct.lattice.abc[0]*10**(-8)
    for i in range(3):
        if not round(struct.lattice.angles[i])==90:
            data['a'] = struct.lattice.abc[(i+1)%3]*10**(-8)
            data['c'] = struct.lattice.abc[i]*10**(-8)
    return data

inp = os.path.join(os.getcwd(),'input.inp')
vdir = get_freq_name(inp,'v')
Edir = get_freq_name(inp,'E')

for i in vdir.keys(): 
    vdir_denom[i] = os.path.join(os.getcwd(),'phonon_'+vdir[i][1]+'_parse')
    vdir_num[i] = os.path.join(os.getcwd(),'phonon_'+vdir[i][0]+'_parse')
for i in Edir.keys():
    Edir_neb[i] = os.path.join(os.getcwd(),'neb_'+Edir[i][1]+'_stat')
    Edir_def[i] = os.path.join(os.getcwd(),'defect_'+Edir[i][0]+'_stat')

def get_HVf(per_dir,def_dir,numatom):
    perfect = open(os.path.join(per_dir,'OSZICAR'),'r')
    defect = open(os.path.join(def_dir,'OSZICAR'),'r')
    for line in perfect:
        if '1 F' in line:
            ene_perf = float(getinfo(line)[2])
    for line in defect:
        if '1 F' in line:
            ene_def = float(getinfo(line)[2])
    perfect.close()
    defect.close()
    HVf = ene_def - (numatom - 1)*ene_perf/numatom
    return HVf

def get_saddle(Edir,Edir_neb):
    for freq in Edir.keys():
        energy = []
        images = 1
        for folder in os.listdir(Edir_neb[freq]):
            if '0'+str(images+2) in folder:
                images = images + 1               
        getsaddle=[os.path.join(Edir_neb[freq],'0'+str(j),'OSZICAR') for j in range(images+2)]
        for j in range(images+2):
            free = open(getsaddle[j],'r')
            for line in free:
                if '1 F' in line:
                    energy.append(float(getinfo(line)[2]))
            free.close()
        enesaddle[freq] = max(energy)
    return enesaddle

def get_end(Edir,Edir_def):
    for freq in Edir.keys():
	fp = open(os.path.join(Edir_def[freq],'OSZICAR'),'r')
        for line in fp:
	    if '1 F' in line:
		eneend[freq] = float(getinfo(line)[2])
	fp.close()
    return eneend

def get_v(vdir,vdir_num,vdir_denom):
    for freq in vdir.keys():
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
    fn.close()
    fd.close()
    return v


k = 11604.52
for freq in Edir.keys():
    enebarr[freq] = get_saddle(Edir,Edir_neb)[freq] - get_end(Edir,Edir_def)[freq]
v = get_v(vdir,vdir_num,vdir_denom)
numatom = get_latt()['No.']
a = get_latt()['a']
c = get_latt()['c']
HVf = get_HVf(per_dir,def_dir,numatom)


if len(vdir)==5:
    D = []
    f0 = 0.7815
    for t in range(20):
        T=10000./(t+1)
        for freq in vdir.keys():
           jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*k/T)
        alpha = jfreq['v4']/jfreq['v0']
        F_num = 10*np.power(alpha,4) + 180.5*np.power(alpha,3) + 927*np.power(alpha,2) + 1341*alpha
        F_denom = 2*np.power(alpha,4) + 40.2*np.power(alpha,3) + 254*np.power(alpha,2) + 597*alpha + 435
        FX = 1-(1.0/7.0)*(F_num/F_denom)
        f2_num = 1+3.5*FX *(jfreq['v3']/jfreq['v1'])
        f2_denom = 1+(jfreq['v2']/jfreq['v1']) + 3.5*FX*(jfreq['v3']/jfreq['v1'])
        f2 = f2_num / f2_denom       
        Vacconc = np.exp(-k*HVf/T)       
        Dself = f0*Vacconc*a**2*jfreq['v0']
        D.append(np.log10(Dself*f2*jfreq['v2']*jfreq['v4']/f0/jfreq['v0']/jfreq['v3']))
    t=np.linspace(0.1,2.0,20)
    plt.plot(t,D,'.')
    plt.show()
    
elif len(vdir)==8:
    Dperp = []
    Dpara = []
    F = 0.736
    lambda1 = a
    lambda2b = a/(3**0.5)
    HB = -0.1 # TBD!
    for t in range(20):
        T=10000./(t+1)
        for freq in vdir.keys():
           jfreq[freq] = v[freq] * np.exp(-enebarr[freq.replace('v','E')]*k/T)
        f1z = (2*jfreq['vap']+7*F*jfreq['vcp'])/(2*jfreq['vap']+7*F*jfreq['vcp']+2*jfreq['vXp'])
        S1b_num = (3**0.5)*jfreq['vX']/jfreq['va']*lambda1+jfreq['vXp']/jfreq['vap']*lambda2b*(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
        S1b_denom = 3 - (2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])*(2+3*jfreq['vbp']/jfreq['vap']+7*F*jfreq['vcp']/jfreq['vap']+2*jfreq['vXp']/jfreq['vap']-1./(2+3*jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+jfreq['vX']/jfreq['va']))
        S1b = S1b_num / S1b_denom
        S2x=((3**0.5)*S1b-jfreq['vX']/jfreq['va']*lambda1)/(2+jfreq['vb']/jfreq['va']+7*F*jfreq['vc']/jfreq['va']+2*jfreq['vX']/jfreq['va'])
        f1b=1+2*S1b/lambda2b
        f2x=1+2*S2x/lambda1
        C2=np.exp(-k*(HVf+HB)/T)
        Dperp.append(np.log10(a**2*(3*f2x*jfreq['vX']+f1b*jfreq['vXp'])*C2/2))
        Dpara.append(np.log10(c**2*f1z*jfreq['vXp']*C2*0.75))
    t=np.linspace(0.1,2.0,20)
    plt.plot(t,Dperp,'.',t,Dpara,'+')
    plt.show()
