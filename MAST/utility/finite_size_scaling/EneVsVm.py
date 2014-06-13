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
from pymatgen.io import vaspio

import GenSC

class writer :
        def __init__(self, *writers) :
                self.writers = writers

        def write(self, text) :
                for w in self.writers :
                        w.write(text)
                        
def psTotEne(oszicar):
    dummyOszcar = vaspio.vasp_output.Oszicar(oszicar).electronic_steps
    TotEne = dummyOszcar[len(dummyOszcar)-1][len(dummyOszcar[len(dummyOszcar)-1])-1]['E']
    return TotEne
    
def CalcV_M(structure):
    dummyPoscarLatStr = mg.core.structure.Structure(structure.lattice,["F-"],[[0,0,0]])
    V_M = ewald.EwaldSummation(dummyPoscarLatStr).total_energy*(-2)
    return V_M


        
def collectEformV_M(perfectDir,defectDirs,defectChg,mu):
    #mu is the chemical potential needed 
    
    Eform = [None]*len(defectDirs)
    V_M = [None]*len(defectDirs)
 
    prf=mg.io.vaspio.vasp_output.Vasprun(perfectDir+"/vasprun.xml")

    prf_ene=prf.final_energy
    #TotEnePr = psTotEne("OSZICAR")

    prf_strct=prf.structures[len(prf.structures)-1]
    
    prf_comps=prf_strct.composition.element_composition
    #perfect_comps = mg.io.vaspio.Poscar.from_file("CONTCAR").structure.composition    
    prf_natom = prf_comps.num_atoms
    
    prf_elmnts = prf_comps.elements
    
    #gap,cbm,vbm,ifdirectgap=prf.eigenvalue_band_properties
    #print (defectDirs[i]+" "+str(V_M[i])+" "+str(gap)+" "+str(cbm)+" "+str(vbm))

    
    all_elmnts = prf_elmnts               
    for i in range(len(defectDirs)):
        
        dfct=mg.io.vaspio.vasp_output.Vasprun(defectDirs[i]+'/vasprun.xml')

        dfct_strct=dfct.structures[len(dfct.structures)-1]
        #defct_strct=mg.io.vaspio.Poscar.from_file("CONTCAR").structure
        
        V_M[i] = CalcV_M(dfct_strct)

        dfct_comps=dfct_strct.composition.element_composition    
        #defct_comps=defct_strct.composition             

        dfct_natom = dfct_comps.num_atoms
        
        dfct_ene = dfct.final_energy
        #defct_ene = psTotEne("OSZICAR")
        
        scalingF = round(dfct_natom/prf_natom)
                    
        dfct_elmnts = dfct_comps.elements
        #print dfct_elmnts
        
        if i == 0:

            for dummy_elmnt1 in dfct_elmnts:
                if dummy_elmnt1 not in all_elmnts:
                    all_elmnts.append(dummy_elmnt1)

            all_amnt_change=[None]*len(all_elmnts) 
            
            for j, dummy_elmnt2 in enumerate(all_elmnts):
                if dummy_elmnt2 in prf_elmnts:
                    dummy_prf_amnt = prf_comps.__getitem__(dummy_elmnt2)
                else:
                    dummy_prf_amnt = 0
                if dummy_elmnt2 in dfct_elmnts:
                    dummy_dfct_amnt = dfct_comps.__getitem__(dummy_elmnt2)
                else:
                    dummy_dfct_amnt = 0
                
                all_amnt_change[j]=scalingF*dummy_prf_amnt-dummy_dfct_amnt                                    

        Eform[i] =  dfct_ene - scalingF*prf_ene 

        #Eform[i] += defectChg*(vbm+mu['e']*gap)
        Eform[i] += defectChg*mu['e']

        for k in range(len(all_elmnts)):

            if all_amnt_change[k] != 0:

                if all_elmnts[k].symbol in mu:                    
                    Eform[i] += mu[all_elmnts[k].symbol]*all_amnt_change[k]
                else:
                    raise RuntimeError("Chemical potential of "+all_elmnts[k].symbol
                                        +"is not provided!")       


    return(V_M,Eform)

def linearFit(listX,listY):
    listXnp = np.array(listX)
    listYnp = np.array(listY)
    slope, intercept, r_value, p_value, std_err = stats.linregress(listXnp,listYnp)
    
    print "The results of linear fitting:"
    print "slope: ", slope
    print "intercept: ", intercept
    print "correlation coefficent r: ", r_value
    print "standard deviation: ", std_err
    print "" 
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
    plt.hold(False) 
    
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
    #fig.text(0.35, 0.25, fittedEq) 
    
    fig.text(0.14,0.25,("Slope: " + str(round(slpIntcpt[0],2))))       
    fig.text(0.14,0.21,("Intercept: " + str(round(slpIntcpt[1],2))))
    fig.text(0.14,0.17,("Correlation coefficent r: " + str(round(slpIntcpt[2],2))))       
    fig.text(0.14,0.13,("Standard deviation: " + str(round(slpIntcpt[3],2)))) 
    
    for i, dirName in enumerate(dirNames):
        fig.text(0.12+0.8*(listX[i]-Left)/(Right-Left)-0.03,
                 0.10+0.8*(listY[i]-Bottom)/(Top-Bottom)+0.01,dirName)       
    
    fig.suptitle(figTitle, fontsize=12)
                          
    plt.show()    

if __name__ == "__main__":
    from MAST.recipe import recipeutility as ru
    fp = open('madelung_utility/V_M.txt','r').readlines()
    plotTitle='Madelung Plot'
    defectDirs=[]
    mu=dict()
    for line in range(len(fp)):
        if fp[line].strip('PlotTitle')[0].strip()=='': plotTitle=fp[line].split('PlotTitle')[1].strip()
        if fp[line].strip('perfect')[0].strip()=='': perfectDir=fp[line].split()[1]  
        if fp[line].strip('mu')[0].strip()=='': 
            labels=fp[line].split('mu')[1].split(',')
            for i in range(len(labels)):
                mu[labels[i].split(':')[0].strip()]=float(labels[i].split(':')[1].strip())
    '''
    fp = open('personal_recipe.txt','r').readlines()
    for line in range(len(fp)):
        if 'madelung_utility' in fp[line]:
            defectDirs.append(fp[line-1].split()[0])
    '''
    defectDirs = ru.read_recipe('personal_recipe.txt')[1]['madelung_utility']
    saved = sys.stdout
    fout = file('madelung_utility/out.log', 'w')
    sys.stdout = writer(sys.stdout, fout)
    
    sc_struct=mg.io.vaspio.Poscar.from_file(defectDirs[0]+'/POSCAR').structure
    sc_potcar=mg.io.vaspio.Potcar.from_file(defectDirs[0]+'/POTCAR')
    sc_incar=mg.io.vaspio.Incar.from_file(defectDirs[0]+'/INCAR') 
    defchg=GenSC.defChg(sc_struct,sc_potcar,sc_incar)
    print "defect charge: ", defchg
            
    (V_M,Eform_orig)=collectEformV_M(perfectDir,defectDirs,defchg,mu)

    slpintcpt=linearFit(V_M,Eform_orig) 
    FS_corr=-np.array(V_M)*slpintcpt[0]
    Eform_corr=np.array(Eform_orig)+FS_corr

    print ("Folder " " V_M " "  EForm_orig" "  EForm_corrct " " FSCorrection")
    print ("-----------------------------------------------------")
    for j in range(len(defectDirs)):
        print (defectDirs[j]+"   "+
               str(round(V_M[j],2))+"   "+
               str(round(Eform_orig[j],5))+"    "+
               str(round(Eform_corr[j],5))+"       "+
               str(round(FS_corr[j],5)))    
    print ""       

    dielectConst=-defchg**2/2/slpintcpt[0]
    print "effective permittivity epsilon_fit: ", round(dielectConst,2)

    plotFit(V_M,Eform_orig,'Madelung Potential v_M (eV)','Formation Energy (eV)',slpintcpt,plotTitle,defectDirs)
    
    sys.stdout = saved
    fout.close()
