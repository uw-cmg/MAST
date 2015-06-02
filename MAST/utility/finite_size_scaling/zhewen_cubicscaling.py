#!/usr/bin/python
import sys, getopt
import numpy as np
from scipy.optimize import curve_fit
import pandas
from MAST.utility.finite_size_scaling.parse_input import ParseInput

class InputError(Exception):
    pass
               
class CubicScaling:
    
    def __init__(self,inputfile,method):
        self.inputfile = inputfile
        self.method = ParseInput().read_input()['method']

    def sortQ(self,Q):
        for i in range(len(Q)):
            if Q[i][-2]=='p': Q[i]=int(Q[i][-1])
            elif Q[i][-2]=='n': Q[i]=-int(Q[i][-1])
        Q.sort()
        for i in range(len(Q)):
            if Q[i]<0: Q[i]='q=n'+str(abs(Q[i])) 
            else: Q[i]='q=p'+str(Q[i])
        return Q
    
    def CubicFit(self,func,listX,listY):
        params = curve_fit(func,listX,listY)[0]
        if func==self.L3:
            SSE = sum((listY-func(listX,params[0],params[1]))**2)
        elif func==self.L1L3:
            SSE = sum((listY-func(listX,params[0],params[1],params[2]))**2)
        SST = sum((listY-np.mean(listY))**2)
        r_coeff = np.sqrt(1 - SSE/SST)
        return [params,r_coeff]
    
    def L3(self,L,a,b):
        return a*L**(-3) + b
    
    def L1L3(self,L,a,b,c):
        return a*L**(-3) + b*L**(-1) + c
    
    def FitStats(self,method,order,L,DFE):
        if len(L)>=order:
            params = self.CubicFit(method,L,DFE[0])
            results = data = "%.4f"%params[0][-1]
        elif len(L)==1:
            results = data = "%.4f"%DFE[0][0]
        else:
            results = data = 'NaN'
        if DFE[-1]==None: r_header = ['E0']
        else: r_header = ['E0_GGA','E0_HSE']
        if DFE[-1]=='NaN': results += " NaN"
        else:
            r_header = ['E0_GGA','E0_HSE']
            if results == 'NaN': results += " NaN"
            else: results += " %.4f"%(float(results)+DFE[-1])
        stats = ''
        warning = ''
        if len(L)>=order:
            if method == self.L3:
                stats = "%.4f  %s  %.4f  "%(params[0][0],data,params[1])
            elif method == self.L1L3:
                stats = "%.4f  %.4f  %s  %.4f  "%(params[0][0],params[0][1],data,params[1])
            if len(L)==order:
                warning += 'ExactFit;'
            if params[1] < 0.9:
                warning += '|R|<0.9;'
            if abs(DFE[0][-1]-params[0][-1]) > 1.0:
                warning += '|E-E0|>1eV;'
            if abs(DFE[0][-1]-params[0][-1])/DFE[0][-1] > 0.2:
                warning += '|E-E0|/E>20%;'
        else:
            stats += data
            warning += ' NoFit;'
        if not warning == '':
            stats += warning
        else: stats += '--'
        if method == self.L3:
            s_header = ['A','E0','|R|','WARNINGS']
        elif method == self.L1L3:
            s_header = ['A1','A2','E0','|R|','WARNINGS']
        
        return [r_header,results,s_header,stats]


    def main(self):
        fp = open(self.inputfile,'r')
        lines = fp.readlines()
        Eform = dict()
        flag = 0
        for i in lines:
            if len(i.split())==1:
                N = i.strip()
                Eform[N] = dict()
            elif ':' in i:
                Q = i.split(':')[0].strip() 
                Eform[N][Q] = dict()
                energy = i.split(':')[1].strip().split()
                for j in range(len(size)):
                    if energy[j]=='NaN':
                        Eform[N][Q][size[j]] = 'NaN'
                    else:
                        Eform[N][Q][size[j]] = float(energy[j])
            elif flag==0:
                size = []
                s = i.split()
                for sz in s:
                    if (not 'delta' in sz) and (not sz.split('x')[0].isdigit()):
                        raise InputError('The input scale size must be in the format like 2x*!')
                    else:
                        size.append(sz)
                flag = 1
        f_results = open('%s_%s.txt'%(self.inputfile.split('_')[0],self.method),'w')
        f_stats = open('%s_%s.stat'%(self.inputfile.split('_')[0],self.method),'w')
        f_stats.write("Statistics for DFE finite size scaling\n")
        if self.method == 'L3':
            f_stats.write("Method: E = A/L^3 + E0\n\n")
        elif self.method == 'L1L3':
            f_stats.write("Method: E = A1/L^3 + A2/L + E0\n\n")
        
        for N in Eform.keys():
            f_results.write("%s\n"%N)
            f_stats.write("%s\n"%N)
            Results = []
            Stats = []
            charge = self.sortQ(Eform[N].keys())
            for Q in charge:
                DFE = [np.array([]),None]
                L = np.array([])
                for S in size:
                    if Eform[N][Q][S] == 'NaN' and (not S=='delta_E'): pass
                    else:
                        if 'delta' in S.lower():
                            DFE[1] = Eform[N][Q][S]
                        else:
                            DFE[0] = np.append(DFE[0],Eform[N][Q][S])
                            L = np.append(L,float(S.split('x')[0]))
                if self.method == 'L3':
                    [r_header,results,s_header,stats] = self.FitStats(self.L3,2,L,DFE)
                elif self.method == 'L1L3':
                    [r_header,results,s_header,stats] = self.FitStats(self.L1L3,3,L,DFE)
                Results.append(results.split())
                Stats.append(stats.split())
            for q in range(len(charge)): charge[q] += ':'
            Data = pandas.DataFrame(Stats,charge,s_header)
            maximum = max([8,Data['WARNINGS'].str.len().max()])
            s_header[-1] = '{{:^{}s}}'.format(maximum).format('WARNINGS')
            f_results.write(pandas.DataFrame(Results,charge,r_header).to_string()+'\n\n')
            f_stats.write(pandas.DataFrame(Stats,charge,s_header).to_string(formatters={s_header[-1]:'{{:<{}s}}'.format(maximum).format})+'\n\n')

if __name__ == "__main__":
    try: opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'cubicscaling.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'cubicscaling.py -i <inputfile>'
            sys.exit()
        elif opt == "-i":
            inputfile = arg
    CubicScaling(inputfile=inputfile,method=ParseInput().read_input()['method']).main()
