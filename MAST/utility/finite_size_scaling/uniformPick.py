#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-05-16
##############################################################
import sys
import math
import warnings
from collections import OrderedDict

class writer :
        def __init__(self, *writers) :
                self.writers = writers

        def write(self, text) :
                for w in self.writers :
                        w.write(text)

def roof_mean(intA,intB):
    """
    This is a helper method that calculates the integer mean of two integers 
    rounded to roof if the direct mean of the two integers are not a integer.
    """
    if (intA+intB)%2 == 0:
        return (intA+intB)/2
    else:
        return (intA+intB+1)/2

def branchingLR(inputList,outputList,m):
    """
    inputList must be a lexicographically ordered list of real numbers with no
    duplicates, and len(inputList)>=3. 
    Note when len(inputlist)==3, the output for left and right branches will be
    the same. 
    m is the number of expected elemented to be picked, m >=3. 
    """    
    if len(inputList)<3:
        raise RuntimeError("The input list only has "+str(len(inputList))+
        " elements, but this algorithm requires it to have >=3 elements.")  
    elif m<3:
        raise RuntimeError("Your input paramter m is "+str(m)+
        " , but this algorithm requires it to be >=3.")     
            
    listMin=inputList[0]
    listMax=inputList[len(inputList)-1] 
           
    Step=(listMax-listMin)/(m-1)

    ####---------min end---------------####              
    for i in range(1,len(inputList),1):
        if inputList[i]>listMin+Step:
            break
            
    if i==1 or abs(inputList[i]-(listMin+Step))<abs(inputList[i-1]-(listMin+Step)):
        leftPick=inputList[i]
        newInputListL=inputList[i:len(inputList)]
    else:        
        leftPick=inputList[i-1]
        newInputListL=inputList[(i-1):len(inputList)]
        
    newOutputListL=list(outputList)
    newOutputListL.append(leftPick)     
        
    ####---------max end---------------####             
    for i in range((len(inputList)-2),0,-1):
        if inputList[i]<listMax-Step:
            break

    if i==len(inputList)-2 or abs(inputList[i]-(listMax-Step))<abs(inputList[i+1]-(listMax-Step)):
        rightPick=inputList[i]
        newInputListR=inputList[0:(i+1)]
    else:
        rightPick=inputList[i+1]
        newInputListR=inputList[0:((i+1)+1)] 

    newOutputListR=list(outputList)
    newOutputListR.append(rightPick)     
        
    return newInputListL,newOutputListL,newInputListR,newOutputListR

def pickSubList(inputList,m):
    
    
    if len(inputList)<m:
        raise RuntimeError("The input list only has "+str(len(inputList))+
        " elements, but you asked to pick "+str(m)+" elements out of the list!")
    elif len(inputList)==m:
        return sorted(list(OrderedDict.fromkeys(inputList))) 
    elif len(inputList)<3:
        raise RuntimeError("The input list only has "+str(len(inputList))+
        " elements, but this algorithm requires it to have >=3 elements.")        
    elif type(m) != int or m<3:
        raise RuntimeError("You asked to select "+str(m)+
        " elements, but this algorithm requires it to be >=3.")        
      

    #removes duplicates from and order the list      
    inputList=sorted(list(OrderedDict.fromkeys(inputList)))
    
    outputList=[inputList[0],inputList[len(inputList)-1]] 

    allInputListDic={}
    allOutputListDic={}
    
    allInputListDic[str(0)+str(0)]=inputList
    allOutputListDic[str(0)+str(0)]=outputList
            
    for i in range(m-2):
        #print "i is "+ str(i)
        for j in range(2**i):
            #print "j is "+ str(j) 
            if len(allInputListDic[str(i)+str(j)])>=3:           
                (newInputListL,newOutputListL,
                 newInputListR,newOutputListR)=branchingLR(allInputListDic[str(i)+str(j)],
                                                       allOutputListDic[str(i)+str(j)],m-i)
            else:
                newInputListL=allInputListDic[str(i)+str(j)]
                newOutputListL=allOutputListDic[str(i)+str(j)]
                newInputListR=allInputListDic[str(i)+str(j)]
                newOutputListR=allOutputListDic[str(i)+str(j)]             
            
            allInputListDic[str(i+1)+str(2*j+0)]=newInputListL
            allOutputListDic[str(i+1)+str(2*j+0)]=newOutputListL
            allInputListDic[str(i+1)+str(2*j+1)]=newInputListR
            allOutputListDic[str(i+1)+str(2*j+1)]=newOutputListR         
  
    bestUniformity=sys.float_info.max
    bestOutputList=[]
    idealSpacing=(inputList[len(inputList)-1]-inputList[0])/(m-1)
    i=m-2
    while bestOutputList==[]:
        for j in range(2**(i-1)):
            for k in [0,1]:
                dummyList=allOutputListDic[str(i)+str(2*j+k)]
                if len(dummyList) == m:
                    dummyUniformity=0
                    for l in range(len(dummyList)-1):
                        dummyUniformity=dummyUniformity+(
                            abs(dummyList[l+1]-dummyList[l])-idealSpacing)**2
                    dummyUniformity=math.sqrt(dummyUniformity/(len(dummyList)-1))                    
                    if dummyUniformity<bestUniformity:
                        bestOutputList=dummyList
        i=i-1
    
    if len(bestOutputList)<m:
        warnings.warn("You asked to pick "+str(m)+" elments out of the inputlist, \
                       but this algorithm only managed to pick"+str(len(bestOutputList))+
                      " elements. The intputlist is probably a very non-uniformly \
                       distributed list.")
                    
    return sorted(bestOutputList)    
 

if __name__ == "__main__":

    saved = sys.stdout
    fout = file('out.log', 'w')
    sys.stdout = writer(sys.stdout, fout)

    testlist=[9.67,6.16,2.59,-0.98,5.0,3.24,1.46,2.24,1.07,0.27,4.84, 
              3.97,3.08,3.48,2.9,2.5,3.22,2.84,2.42]
    
    for n in range(3,21):
        print pickSubList(testlist,n)
         
    sys.stdout = saved
    fout.close()   
