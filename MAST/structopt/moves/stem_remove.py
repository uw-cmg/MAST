import random
from MAST.structopt.tools.StemCalc import ConvStem
from ase import Atoms,Atom
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 

def stem_remove(indiv, Optimizer):
    atms = indiv[0].copy()
    pixelshift = [0,0,0]
    if Optimizer.stemcalc.parameters['Pixelshift'] == True:
       points = [-0.5993457/2,0,0.5993457/2]
       chisq = numpy.zeros([3,3],dtype=float)
       for x in range(3):
          for y in range(3):
             pixelshift = [points[x],points[y],0]
             simfun=Optimizer.stemcalc.get_image(Optimizer.stemcalc.psf,atms,Optimizer.stemcalc.parameters['Slice size'],Optimizer.stemcalc.parameters['Pixels'],pixelshift)
             chisq[x][y]=Optimizer.stemcalc.compare_functions(Optimizer.stemcalc.expfun,simfun)
       i,j = numpy.unravel_index(chisq.argmin(),chisq.shape)
       pixelshift = [points[i],points[j],0]
    simfun = Optimizer.stemcalc.get_image(Optimizer.stemcalc.psf,atms,Optimizer.stemcalc.parameters['Slice size'],Optimizer.stemcalc.parameters['Pixels'], pixelshift) 
    expfun = Optimizer.stemcalc.expfun

    xmax = Optimizer.stemcalc.parameters['Slice size']
    com = atms.get_center_of_mass()
    #com = [ 44.40963074 , 44.65497562 , 44.90406073]
    #com = numpy.array(com)
    #com += [-0.149836425, 0.29967285, 0]  #pixelshift
    cop = xmax/2.0 
    trans = [cop-i for i in com]
    atms.translate(trans)
    R = atms.arrays['positions']
    ax=[]
    ay=[]
    az=[] 
    for o,t,h in R:
        ax.append(o)
        ay.append(t)
        az.append(h)
    ax = numpy.array(ax)
    ay = numpy.array(ay)
    az = numpy.array(az)

    nx = Optimizer.stemcalc.parameters['Pixels']
    ny = nx
    dx = xmax/nx
    dy = dx
    ix = numpy.array([math.floor(axi/dx) for axi in ax])
    iax = numpy.array([int(math.fmod(iaxi,nx)) for iaxi in ix])
    ibx = numpy.array([int(math.fmod(iaxi+1,nx)) for iaxi in ix])
    iy = numpy.array([math.floor(ayi/dy) for ayi in ay])
    iay = numpy.array([int(math.fmod(iayi,ny)) for iayi in iy])
    iby = numpy.array([int(math.fmod(iayi+1,ny)) for iayi in iy])

 
# construct a map from atom to pixel grid   
    map_atom_pixel = []

    for i in range(nx) :
       map_atom_pixel.append([])
       for j in range(ny) :
          map_atom_pixel[i].append([])

    for iatom in range(len(ix)):
       for x in [iax[iatom],ibx[iatom]] :
          for y in [iay[iatom],iby[iatom]] :
              map_atom_pixel[x][y].append(iatom)
 
#remove atoms from postive delta STEM intensity 
    delta = simfun-expfun
    for i in range(len(expfun)):
        for j in range(len(expfun[0])):
           if expfun[i][j] <= -1000 : #invalid exp points
             delta[i][j] = 0 
    remove_tol = delta.max()*0.3
    maxmove = 1

# construct the pixel list where simfun differs from expfun
    remove_list_x = []
    remove_list_y = []

    for x in range(nx) :
       for y in range(ny) :
          if delta[x][y] > remove_tol: # and x < 40 and x >10 and y <120 and y >60:
             remove_list_x.append(x)
             remove_list_y.append(y)
    remove_list_x = numpy.array(remove_list_x)
    remove_list_y = numpy.array(remove_list_y)

 
#removing atom from delta positive col 
    removed = []
    for imove in range(maxmove): 
       removecol=random.randint(0,len(remove_list_x)-1)
       x_pixel = remove_list_x[removecol]
       y_pixel = remove_list_y[removecol]
       #print 'remove',x_pixel,y_pixel,delta[x_pixel][y_pixel],map_atom_pixel[x_pixel][y_pixel]
       if len(map_atom_pixel[x_pixel][y_pixel]) > 0:
         ind = map_atom_pixel[x_pixel][y_pixel][0]
         atomind_max = ind
         atomind_min = ind
         Zmax = R[ind][2]-cop
         Zmin = R[ind][2]-cop
         for j in range(1,len(map_atom_pixel[x_pixel][y_pixel])):
           ind = map_atom_pixel[x_pixel][y_pixel][j]
           if R[ind][2]-cop > Zmax :
             atomind_max = ind
             Zmax = R[ind][2]-cop
           elif R[ind][2]-cop < Zmin : 
             atomind_min = ind
             Zmin = R[ind][2]-cop  
         atomind_remove = random.choice([atomind_max,atomind_min])      
         removed.append(atomind_remove)
         #print 'remove',atomind_remove,R[atomind_remove][2]
       # remove this atom index from removecol of map_atom_pixel  
         for x in [iax[atomind_remove],ibx[atomind_remove]] :
            for y in [iay[atomind_remove],iby[atomind_remove]] :
               map_atom_pixel[x][y].remove(atomind_remove)
    removed = list(set(removed))
    removed.sort(reverse=True)
    for i in range(len(removed)):
        indiv[0].pop(removed[i]) 

    muttype='STEM'+repr(maxmove)
    Optimizer.output.write('STEM remove performed on individual\n')
    Optimizer.output.write('Index = '+repr(removed)+'\n')
    #Optimizer.output.write(repr(indiv[0])+'\n')
    #Optimizer.output.write(repr(indiv[0].get_positions())+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype

    return indiv
