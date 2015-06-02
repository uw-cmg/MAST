import random
from MAST.structopt.tools.StemCalc import ConvStem
from ase import Atoms, Atom
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 
from collections import Counter

def stem_shift(indiv, Optimizer):
    """Shift n cols of atoms upward/downward together -- this would only change E but not STEM intensity
    Input:
        indiv = individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Output:
        indiv = Altered Individual class object
    """

    atms = indiv[0].copy()

    xmax = Optimizer.stemcalc.parameters['Slice size']
    atm_num = atms.get_number_of_atoms()
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
    #print 'stem_shift'
    maxshift = 20
    shift_refine = 4

# construct the list of atom to be shifted
    shiftatmlist = []
    flag_shiftcol = False
    while flag_shiftcol == False :
       shiftcol_x = random.randint(0,nx-1)
       shiftcol_y = random.randint(0,ny-1)

       if len(map_atom_pixel[shiftcol_x][shiftcol_y]) > 5 :
           flag_shiftcol = True

    #print 'shiftcol_x,shiftcol_y',shiftcol_x,shiftcol_y

    for dx in range(-maxshift/2,maxshift/2+1) :
      for dy in range(-maxshift/2,maxshift/2+1) :
         if dx**2 + dy**2 <= (maxshift/2)**2 : 
            shiftatmlist.extend(map_atom_pixel[shiftcol_x+dx][shiftcol_y+dy])
    #refinement: to include atoms on circle boundary, avoid choosing part of col.
    for dx in range(-(maxshift+shift_refine)/2,(maxshift+shift_refine)/2+1) :
      for dy in range(-(maxshift+shift_refine)/2,(maxshift+shift_refine)/2+1) :
          if dx**2 + dy**2 > (maxshift/2)**2 and dx**2 + dy**2 <= ((maxshift+shift_refine)/2)**2 :    
            if len(list((Counter(map_atom_pixel[shiftcol_x+dx][shiftcol_y+dy]) & Counter(shiftatmlist)).elements())) > 0 :
               shiftatmlist.extend(map_atom_pixel[shiftcol_x+dx][shiftcol_y+dy])

    shiftatmlist = sorted(set(shiftatmlist))
    #print 'indiv.history_index',indiv.history_index
    #print 'shiftatmlist'
    #for iatom in shiftatmlist:
        #print iatom,numpy.array(R[iatom])

    sign = random.choice([1,-1])
    #sign = -1
    for iatom in shiftatmlist :
       R[iatom][2] += sign * 2.9
     #  print iatom, numpy.array(R[iatom])

    indiv[0] = atms
    muttype='AtmShift'+repr(maxshift)
    Optimizer.output.write('AtmShift performed on individual\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype

    return indiv

