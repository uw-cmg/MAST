#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-06-13 by Zhewen Song
##############################################################
import sys, getopt, os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io.vasp import Oszicar, Outcar, Poscar

import GenSC
'''
class writer :
        def __init__(self, *writers) :
                self.writers = writers

        def write(self, text) :
                for w in self.writers :
                        w.write(text)
'''                        
def psTotEne(oszicar):
    dummyOszcar = Oszicar(oszicar).electronic_steps
    TotEne = dummyOszcar[len(dummyOszcar)-1][len(dummyOszcar[len(dummyOszcar)-1])-1]['E']
    return TotEne
    
def CalcV_M(structure):
    dummyPoscarLatStr = mg.core.structure.Structure(structure.lattice,["F-"],[[0,0,0]])
    V_M = ewald.EwaldSummation(dummyPoscarLatStr).total_energy*(-2)
    return V_M

def linearFit(listX,listY):
    listXnp = np.array(listX)
    listYnp = np.array(listY)
    slope, intercept, r_value, p_value, std_err = stats.linregress(listXnp,listYnp)
    return [slope,intercept,r_value,std_err]   
    
def plotFit(listX,listY,Xlabel,Ylabel,slpIntcpt,figTitle,dirNames):
    listXnp = np.array(listX)
    listYnp = np.array(listY)
    listYnpCorrectted = listYnp-slpIntcpt[0]*listXnp

    Left=min(listX)-0.5
    Right=max(listX)+0.5
        
    listXnpFitted=np.arange(Left,Right,(Right-Left)/6)
    listXnpFitted=np.append(listXnpFitted,Right)
    listYnpFitted = slpIntcpt[0]*listXnpFitted + slpIntcpt[1]
    listYnpCorrecttedFitted = listXnpFitted*0+slpIntcpt[1]
    
    fig = plt.figure()
    ax = fig.add_axes([0.12,0.1,0.8,0.8])
    
    #plt.plot(listXnp,listYnp,'ro',listXnpFitted,listYnpFitted,'r-',
    #         listXnp,listYnpCorrectted,'bo',listXnpFitted,listYnpCorrecttedFitted,'b--')
    line1=plt.plot(listXnp,listYnp,'ro',label='original')
    plt.hold(True)
    line2=plt.plot(listXnpFitted,listYnpFitted,'r-')
    line3=plt.plot(listXnp,listYnpCorrectted,'bo',label='corrected')
    line4=plt.plot(listXnpFitted,listYnpCorrecttedFitted,'b--')
    listXnpFitted,listYnpFitted,'r-',
    #plt.hold(False) 
    
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)

    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.xlim(Left,Right)
    
    dummyBottom=min(min(listYnp),min(listYnpCorrectted))
    dummyTop=max(max(listYnp),max(listYnpCorrectted))
    Bottom=dummyBottom-0.2*(dummyTop-dummyBottom)
    Top=dummyTop+0.2*(dummyTop-dummyBottom)
    plt.ylim(Bottom,Top)

    plt.axvline(x=0,color='b',ls='dashed')
    
    #fittedEq = (Ylabel + " = " + str(slope) + " * " + Xlabel + " + " + str(intercept))
    #fittedEq = ("Y = (" + str(slpIntcpt[0]) + ") * X" + " + (" + str(slpIntcpt[1]) + ")")
    #fittedEq = ("Slope is " + str(round(slpIntcpt[0],2)) + 
    #            ". Intercept is " + str(round(slpIntcpt[1],2))+".")
    
    fig.text(0.14,0.25,("Slope: " + str(round(slpIntcpt[0],2))))       
    fig.text(0.14,0.21,("Intercept: " + str(round(slpIntcpt[1],2))))
    fig.text(0.14,0.17,("Correlation coefficent R: " + str(round(slpIntcpt[2],2))))       
    fig.text(0.14,0.13,("Standard deviation: " + str(round(slpIntcpt[3],2)))) 
    
    # for i, dirName in enumerate(dirNames):
    #    fig.text(0.12+0.8*(listX[i]-Left)/(Right-Left)-0.03,
    #             0.10+0.8*(listY[i]-Bottom)/(Top-Bottom)+0.01,dirName)       
    
    fig.suptitle(figTitle, fontsize=12)
                          

if __name__ == "__main__":
    from MAST.parsers import InputParser
    from MAST.recipe import recipeutility as ru
    from MAST.utility.defect_formation_energy import DefectFormationEnergy as DFE
    import shutil, os
  

    print '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '% Defect Formation Energy vs. Madelung Potential %'
    print '%    from MAST (MAterials Simulation Toolbox)    %'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
 
    dfe = DFE('..')
    fp = open('DFE_vs_VM.log', 'w')
    v_m = dict()
    Eform = dict()
    scalingsize = dfe.input_options.get_item('structure','scaling')
#    ingredients = dfe.recipe_plan.ingredients.keys()
    for size_label in scalingsize.keys():
        dfe = DFE('..')
        dfe._calculate_defect_formation_energies("[%s]"%scalingsize[size_label][0])
        Ef = dfe.get_defect_formation_energies("[%s]"%scalingsize[size_label][0])
        Eform[size_label] = Ef
        for folder in ingredients:
            if (size_label in folder) and ('perfect' in folder):
                struct = Poscar.from_file('../%s/POSCAR'%folder).structure
                v_m[size_label] = CalcV_M(struct)
                break
    sizes = Eform.keys()
    sizes.sort()
    chempot = dfe.input_options.get_item('chemical_potentials')
    conditions = chempot.keys()
    defects = Eform[sizes[0]][conditions[0]].keys()
    defects.sort()
    
    print Eform
    charge_defects = dict()
    for DEF in defects:
        try: charge_defects[DEF.split('_q=')[0]].append(DEF.split('_q='))
        except KeyError: charge_defects[DEF.split('_q=')[0]] = [DEF.split('_q=')]
    efermi=np.linspace(0.0,6.0,num=100)
    for CON in conditions:
        plt.figure()
        for DEF in charge_defects.keys():
            eforms=np.zeros(100)
            for i in range(100):
                intercept = []
                for CHG in charge_defects[DEF]:
                    intercept.append(linearFit(V_M,eform['%s_q=%s'%(DEF,CHG)][CON])+(efermi[i]*float(CHG)))
                eforms[i] = min(intercept)
            plt.plot(efermi,eforms,'-',label='%s'%DEF)
            plt.hold(True)
        #plt.legend()
        plt.xlabel('Fermi Energy eV')
        plt.ylabel('Defect Formation Energy eV')
        plt.savefig('EfvsFermi_%s.png'%(CON))



