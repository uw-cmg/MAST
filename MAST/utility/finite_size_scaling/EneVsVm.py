#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-09-13 by Zhewen Song
##############################################################
import sys, getopt, os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio

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
    dummyOszcar = vaspio.vasp_output.Oszicar(oszicar).electronic_steps
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

    Left=min(listX)-0.2
    Right=max(listX)+0.2
        
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

    plt.legend()

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
    fig.text(0.14,0.17,("Correlation coefficent r: " + str(round(slpIntcpt[2],2))))       
    fig.text(0.14,0.13,("Standard deviation: " + str(round(slpIntcpt[3],2)))) 
    
    # for i, dirName in enumerate(dirNames):
    #    fig.text(0.12+0.8*(listX[i]-Left)/(Right-Left)-0.03,
    #             0.10+0.8*(listY[i]-Bottom)/(Top-Bottom)+0.01,dirName)       
    
    fig.suptitle(figTitle, fontsize=12)
                          

if __name__ == "__main__":
    from MAST.parsers import InputParser
    from MAST.recipe import recipeutility as ru
    from MAST.utility.defect_formation_energy import DefectFormationEnergy as DFE
    import shutil
    
    dfe = DFE('..')
    fp = open('madelung_utility.log', 'w')
    print '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '% Defect Formation Energy vs. Madelung Potential %'
    print '%    from MAST (MAterials Simulation Toolbox)    %'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
    v_m = dict()
    Eform_orig = dict()
    scalingsize = dfe.input_options.get_item('structure','scaling')
    ingredients = dfe.recipe_plan.ingredients.keys()
    for size_label in scalingsize.keys():
        dfe = DFE('..')
        dfe._calculate_defect_formation_energies("[%s]"%scalingsize[size_label][0])
        Ef = dfe.get_defect_formation_energies("[%s]"%scalingsize[size_label][0])
        Eform_orig[size_label] = Ef
        for folder in ingredients:
            if (size_label in folder) and ('perfect' in folder):
                struct = mg.io.vaspio.Poscar.from_file('../'+folder+"/POSCAR").structure
                v_m[size_label] = CalcV_M(struct)
                break
        for folder in ingredients:
            if (size_label in folder) and folder.startswith('defect'):
                shutil.copy('../%s/CONTCAR'%folder,'%s_structure'%folder)
                shutil.copy('../%s/OSZICAR'%folder,'%s_energy'%folder) 
    sizes = Eform_orig.keys()
    sizes.sort()
    chempot = dfe.input_options.get_item('chemical_potentials')
    conditions = chempot.keys()
    defects = Eform_orig[sizes[0]][conditions[0]].keys()
    elements = chempot[conditions[0]]
    string=''
    for ele in elements.keys():
        string+='{0:^8}'.format('mu_%s'%ele)
    for DEF in defects:
        fp.write( "Table for defect %s:\n"%DEF )
        fp.write( "{0:^8}{1:^8}{2:^8}{3:^8}{4:^8}{5:^8}\n".format(' ',"size","V_M","E_orig","FS_corr","E_corr"))
        fp.write(''.join([i*48 for i in '='])+"\n")       
        #Eform = []
        slpintcpt = dict()
        for CON in conditions:
            eform = []
            V_M = []
            for SIZE in sizes:
                eform.append(Eform_orig[SIZE][CON][DEF])
                V_M.append(v_m[SIZE])
            #Eform.append(eform)               
            slpintcpt[CON]=linearFit(V_M,eform) 
            FS_corr=-np.array(V_M)*slpintcpt[CON][0]
            Eform_corr=np.array(eform)+FS_corr
            cases = 'Defect:%s with Conditions:%s'%(DEF,CON)
            plotFit(V_M,eform,'Madelung Potential V_M (eV)','Formation Energy (eV)',slpintcpt[CON],cases,"[%s]"%sizes)    
            plt.savefig("%s_%s.png"%(DEF,CON))
            for i in range(len(V_M)):
                if i==0: blank = CON
                else: blank = " "
                fp.write("{0:^8}{1:^8}{2:^8}{3:^8}{4:^8}{5:^8}\n".format(blank,sizes[i],"%.2f"%V_M[i],"%.2f"%eform[i],"%.2f"%FS_corr[i],"%.2f"%Eform_corr[i]) )
                if i==len(V_M)-1: fp.write(''.join([i*48 for i in '-'])+"\n") 
        fp.write( "{0:^8}{1:^8}{2:^8}{3:^8}{4:^8}{5:^8}\n".format(' ',string,'DFE','slope','R','stdev'))
        fp.write(''.join([i*8*(5+len(elements)) for i in '='])+"\n") 
        for con in conditions:
            mu = ""
            for ele in elements.keys():
                mu+="{0:^8}".format('%.2f'%chempot[con][ele])
            fp.write("{0:^8}{1:^8}{2:^8}{3:^8}{4:^8}{5:^8}\n".format(con,mu,"%.2f"%slpintcpt[con][1],"%.2f"%slpintcpt[con][0],"%.3f"%slpintcpt[con][2],"%.3f"%slpintcpt[con][3]))
            fp.write(''.join([i*8*(5+len(elements)) for i in '-'])+"\n")
        fp.write("\n")
    fp.close()
    print (open('madelung_utility.log','r').read())
