##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-04-25
##############################################################
import sys, getopt, os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio


def psTotEne(oszicar):
    dummyOszcar = vaspio.vasp_output.Oszicar(oszicar).electronic_steps
    TotEne = dummyOszcar[len(dummyOszcar)-1][len(dummyOszcar[len(dummyOszcar)-1])-1]['E']
    return TotEne
    
def CalcV_M(structure):
    dummyPoscarLatStr = mg.core.structure.Structure(structure.lattice,["F-"],[[0,0,0]])
    V_M = ewald.EwaldSummation(dummyPoscarLatStr).total_energy*(-2)
    return V_M
    
def FitnPlot(listX,listY,Xlabel,Ylabel,plotBoolean):
    listXnp = np.array(listX)
    listYnp = np.array(listY)
    slope, intercept, r_value, p_value, std_err = stats.linregress(listXnp,listYnp)

    print 'slope: ', slope
    print 'intercept: ', intercept
    print 'correlation coefficent r: ', r_value
    print 'standard deviation: ', std_err
    
    if plotBoolean:
        listYfitted = slope*listXnp + intercept
        fig = plt.figure()
        ax = fig.add_axes([0.12,0.1,0.8,0.8])
        plt.plot(listXnp,listYfitted,'r-',listXnp,listYnp,'o')
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)

        #fittedEq = (Ylabel + " = " + str(slope) + " * " + Xlabel + " + " + str(intercept))        
        fittedEq = ("Y = (" + str(slope) + ") * X" + " + (" + str(intercept) + ")")
        #fig.text(0.35, 0.75, str('intercept is'+' '+str(intercept)))
        fig.text(0.35, 0.75, fittedEq)
        
        plt.show()    
        

if __name__ == "__main__":
    perfectDir = "perfect_30atoms"
    Dirs = ["2x2x1","2x2x2","2x2x3","3x3x1","4x4x1","3x3x2"]

    FormEneDef = [None]*len(Dirs)
    V_M = [None]*len(Dirs)

    cwdir = os.getcwd()

    os.chdir(perfectDir)
    TotEnePr = psTotEne("OSZICAR")
    NatomPr = mg.io.vaspio.Poscar.from_file("POSCAR").structure.composition.num_atoms
    os.chdir(cwdir)

    for i in range(len(Dirs)):
        os.chdir(Dirs[i])

        dummyTotEneSc = psTotEne("OSZICAR")
        dummmyNatomSc = mg.io.vaspio.Poscar.from_file("POSCAR").structure.composition.num_atoms
        dummyScalF = round(dummmyNatomSc/NatomPr)
        #print dummyScalF
        FormEneDef[i] =  dummyTotEneSc - dummyScalF*TotEnePr
        
        #V_M[i] = CalcV_M("CONTCAR")
        V_M[i] = CalcV_M(mg.io.vaspio.Poscar.from_file("CONTCAR").structure)      
    
        print (Dirs[i]+" "+str(V_M[i])+" "+str(FormEneDef[i]))
    
        os.chdir(cwdir)

    print ""
    FitnPlot(V_M,FormEneDef,'Madelung Potential V_M (eV)','Formation Energy (eV)',True)
    


