import random
from MAST.structopt_stem.tools.StemCalc import ConvStem
from ase import Atoms
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 

def stem_mutation(indiv, Optimizer):
    atms = indiv[0].copy()
    simfun = Optimizer.stemcalc.get_image(Optimizer.stemcalc.psf,atms,Optimizer.stemcalc.parameters['Slice size'],Optimizer.stemcalc.parameters['Pixels'])
    expfun = Optimizer.stemcalc.expfun

    xmax = Optimizer.stemcalc.parameters['Slice size']
    com = atms.get_center_of_mass()
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
 
#move atoms from high STEM intensity to low STEM intensity
    delta = simfun-expfun
    for i in range(len(expfun)):
        for j in range(len(expfun[0])):
           if expfun[i][j] <= -1000 : #invalid exp points
             delta[i][j] = 0 
    remove_tol = delta.max()*0.5
    add_tol = delta.min()*0.5
    maxmove = 5

# construct the pixel list where simfun differs from expfun
    remove_list_x = []
    remove_list_y = []
    add_list_x = []
    add_list_y = []

    for x in range(nx) :
       for y in range(ny) :
          if delta[x][y] > remove_tol:
             remove_list_x.append(x)
             remove_list_y.append(y)
          elif delta[x][y] < add_tol:
             add_list_x.append(x)
             add_list_y.append(y)
    remove_list_x = numpy.array(remove_list_x)
    remove_list_y = numpy.array(remove_list_y)
    add_list_x = numpy.array(add_list_x)
    add_list_y = numpy.array(add_list_y)

 
#loop for moving atom from delta positive col to delta negative col
    moved = []
    for imove in range(maxmove): 
       removecol_atomnum = 0
       while removecol_atomnum < 1:
           removecol=random.randint(0,len(remove_list_x)-1)
           x_pixel = remove_list_x[removecol]
           y_pixel = remove_list_y[removecol]
           removecol_atomnum = len(map_atom_pixel[x_pixel][y_pixel]) 
       print x_pixel, y_pixel, map_atom_pixel[x_pixel][y_pixel]
       ind = map_atom_pixel[x_pixel][y_pixel][0]
       atomind_remove = ind
       Zmax = abs(R[ind][2]-cop)
       for j in range(1,len(map_atom_pixel[x_pixel][y_pixel])):
          ind = map_atom_pixel[x_pixel][y_pixel][j]
          if abs(R[ind][2]-cop) > Zmax :
             atomind_remove = ind
             Zmax = abs(R[ind][2]-cop)
       moved.append(atomind_remove)
       # remove this atom index from removecol of map_atom_pixel  
       map_atom_pixel[x_pixel][y_pixel].remove(atomind_remove)

       # check to see if this negative col contains any atoms already, 
       # if not, skip this col, as it is hard to assign atomic coordinates.  need better algorithm here
       addcol_atomnum = 0
       while addcol_atomnum < 1 :
           addcol = random.randint(0,len(add_list_x)-1)
           x_pixel = add_list_x[addcol]
           y_pixel = add_list_y[addcol]
           addcol_atomnum = len(map_atom_pixel[x_pixel][y_pixel]) 
       if addcol_atomnum > 0 :
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
         atomind_add = random.choice([atomind_max,atomind_min])
         if atomind_add == atomind_max :
            sign = 1
         elif atomind_add == atomind_min :
            sign = -1
       # now add the removed atom index back to addcol of map_atom_pixel 
         map_atom_pixel[x_pixel][y_pixel].append(atomind_remove)        

       # update atomic coord 
       R[atomind_remove][0] = R[atomind_add][0]
       R[atomind_remove][1] = R[atomind_add][1]
       R[atomind_remove][2] = R[atomind_add][2] + sign*2.6

    indiv[0] = atms
    muttype='STEM'+repr(maxmove)
    Optimizer.output.write('STEM Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(moved)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    #Optimizer.output.write(repr(indiv[0].get_positions())+'\n')
    #Optimizer.output.write('moved atoms='+repr(poplist)+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype

    return indiv
