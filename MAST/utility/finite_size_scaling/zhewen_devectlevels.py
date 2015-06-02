import sys, getopt
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import MAST.utility.finite_size_scaling.colors as colors
from MAST.utility.finite_size_scaling.parse_input import ParseInput
import pandas
            
class DefectLevels:
    def __init__(self,inputfile,XC):
        self.inputfile = inputfile
        self.XC =XC
    
    def Q2charge(self,Q):
        if Q.split('=')[1][0] == 'n' or Q.split('=')[1][0] == '-':
            charge = -float(Q.split('=')[1][1:])
        else: charge = float(Q.split('=')[1][1:])
        return charge

    def get_cross(self,chg,i,DFE,N,id):
        k1,k2 = self.Q2charge(chg[i]),self.Q2charge(chg[i+1])
        b1,b2 = DFE[N][chg[i]][id],DFE[N][chg[i+1]][id]
        return [(b2-b1)/(k1-k2),(k1*b2-b1*k2)/(k1-k2)]

    def get_eforms(self,DFE,gap,id):
        bins = 1000
        efermi=np.linspace(0.0,gap,bins)
        eforms=dict()
        cross = dict()
        for N in DFE.keys():
            eforms[N] = np.zeros(bins)
            chg = [None]*1001
            cross[N] = dict()      
            for i in range(bins):
                tmp = dict()
                for Q in DFE[N].keys():
                    if np.isnan(DFE[N][Q][id]): continue
                    charge = self.Q2charge(Q)
                    tmp[Q] = DFE[N][Q][id]+efermi[i]*charge
                chg[i+1] = min(tmp,key=tmp.get)
                eforms[N][i] = tmp[chg[i+1]]
                if i==0:
                    cross[N][chg[i+1]] = [[efermi[i],eforms[N][i]]]
                elif (not chg[i+1]==chg[i]) and i<bins-1:
                    cross[N][chg[i+1]] = [self.get_cross(chg,i,DFE,N,id)]
                    cross[N][chg[i]].append(self.get_cross(chg,i,DFE,N,id))
                elif i==bins-1:
                    cross[N][chg[i+1]].append([efermi[i],eforms[N][i]])
        return [efermi,eforms,cross]
            
    def plot_style(self,defects,efermi,eforms,name):
        plt.figure()
        COLORS = colors.COLORS
        num = 0
        ax = plt.subplot(111)
        for N in defects:
            ax.plot(efermi,eforms[N],color=COLORS[num%len(defects)],linewidth=2.0,label='%s'%N)
            num += 1
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.tick_params(labelsize=18)
        plt.xlim(0,efermi[-1])
        plt.xlabel('Fermi Energy (eV)',size=16)
        plt.ylabel('Defect Formation Energy (eV)',size=16)
        plt.savefig(name)     
        return
    
    def write_txt(self,defects,efermi,eforms,cross,name):
        fw = open(name,'w')
        for e in range(len(efermi)): efermi[e] = '%.3f'%efermi[e]
        Data = pandas.DataFrame(eforms,efermi,defects)
        fw.write(Data.to_string(float_format=lambda x: '%.4f' % x)+'\n\n')
        fw.write('Transition levels\n')
        for N in cross.keys():
            fw.write('%s\n'%N)
            charge = sorted(cross[N],key=cross[N].get)
            ranges = []
            for Q in charge:
                string = '%.4f %.4f %.4f %.4f'%(cross[N][Q][0][0],cross[N][Q][1][0],cross[N][Q][0][1],cross[N][Q][1][1])
                ranges.append(string.split())
            for q in range(len(charge)):    
                charge[q]+=':'
            Data = pandas.DataFrame(ranges,charge,['Ef1','Ef2','DFE1','DFE2'])
            fw.write(Data.to_string()+'\n\n')
        return
                                           
    def stretch(self,efermi,eforms,cross,gap):
        for i in range(len(efermi)):
            efermi[i] = efermi[i]*gap[1]/gap[0]
        for N in cross.keys():
            for Q in cross[N].keys():
                cross[N][Q][0][0] = cross[N][Q][0][0]*gap[1]/gap[0]
                cross[N][Q][1][0] = cross[N][Q][1][0]*gap[1]/gap[0]      
        return [efermi,eforms,cross]
        
    def main(self):
        fp = open(self.inputfile,'r')
        lines = fp.readlines()
        cross = []
        DFE = defaultdict(lambda:defaultdict(defaultdict)) 
        for l in lines:
            if len(l.split())==1:
                N = l.split()[0].strip()
            elif len(l.split())>1 and ':' in l:
                DFE[N][l.split(':')[0].strip()] = np.array(l.split(':')[1].split(),dtype='float')
        for id in range(len(self.XC)):
            [efermi,eforms,crossid] = self.get_eforms(DFE,self.XC[id]['gap']['mygap'][0],id)
            cross.append(crossid)
            name = "%s_%s"%(self.inputfile.split('.')[0],self.XC[id]['tag'])
            self.plot_style(DFE.keys(),efermi,eforms,name+'.eps')
            self.write_txt(DFE.keys(),efermi,eforms,cross[id],name+'.txt')
            if len(self.XC[id]['gap']['mygap'])>=2:
                for j in range(1,len(self.XC[id]['gap']['mygap'])):
                    [efermi,eforms,cross[id]] = self.stretch(efermi,eforms,cross[id],self.XC[id]['gap']['mygap'])
                    name = "%s_%s_stretched%s"%(self.inputfile.split('.')[0],self.XC[id]['tag'],j)
                    self.plot_style(DFE.keys(),efermi,eforms,name+'.eps')
                    self.write_txt(DFE.keys(),efermi,eforms,cross[id],name+'.txt')
                    

if __name__ == "__main__":
    try: opts, args = getopt.getopt(sys.argv[1:],"hi:",["ifile="])
    except getopt.GetoptError:
        print 'defectlevels.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'defectlevels.py -i <inputfile>'
            sys.exit()
        elif opt == "-i":
            inputfile = arg
    inp = ParseInput().read_input()
    XC = [inp['GGA'],inp['HSE']]
    DefectLevels(inputfile=inputfile,XC=XC).main()
