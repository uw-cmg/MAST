import os
try:
    from ase import Atom, Atoms
    from ase.calculators.vasp import Vasp
    from ase.optimize import BFGS
    from ase.units import GPa
    from ase.calculators.neighborlist import NeighborList
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt eval_energy.py, ASE must be installed."
from MAST.structopt_stem.inp_out.write_xyz import write_xyz
from MAST.structopt_stem.tools.setup_calculator import setup_calculator
from MAST.structopt_stem.tools.find_defects import find_defects
from MAST.structopt_stem.tools.check_cell_type import check_cell_type
from MAST.structopt_stem.fingerprinting import get_fingerprint
from MAST.structopt_stem.tools.lammps import LAMMPS
import numpy
import math
try:
    from mpi4py import MPI
except ImportError:
    print "NOTE: mpi4py is not installed. To use certain features in Structopt eval_energy.py, mpi4py must be installed."
import logging
import pdb
import shutil
import time
import scipy
import random

def eval_energy(Optimizer, individ):
    """Function to evaluate energy of an individual
    Inputs:
        input = [Optimizer class object with parameters, Individual class structure to be evaluated]
    Outputs:
        energy, bul, individ, signal
        energy = energy of Individual evaluated
        bul = bulk structure of Individual if simulation structure is Defect
        individ = Individual class structure evaluated
        signal = string of information about evaluation
    """
    #logger = initialize_logger(Optimizer.loggername)
    logger = logging.getLogger(Optimizer.loggername)
    if 'MAST' in Optimizer.calc_method:
        energy = individ.energy
        bul = individ.bulki
        signal = 'Received MAST structure\n'
        logger.info('Received individual index = {0} from MAST with energy {1}. Returning with no evaluation'.format(
            individ.index, individ.energy))
    else:
        if Optimizer.parallel: 
            rank = MPI.COMM_WORLD.Get_rank()
        logger.info('Received individual HI = {0} with energy {1} for energy evaluation'.format(
            individ.history_index, individ.energy))
        STR='----Individual ' + str(individ.history_index)+ ' Optimization----\n'
        indiv=individ[0]
        if 'EE' in Optimizer.debug:
            debug = True
        else:
            debug = False
        if debug: 
            write_xyz(Optimizer.debugfile,indiv,'Received by eval_energy')
            Optimizer.debugfile.flush()
            logger.debug('Writing recieved individual to debug file')
        # Establish individual structure for evaluation.  Piece together regions when necessary.
        if Optimizer.structure=='Defect':
            indi=indiv.copy()
            bulk=individ.bulki
            nat=indi.get_number_of_atoms()
            if debug: 
                logger.info('Extending defect structure to include bulk len(r1+r2)={0} len(bulk)={1}'.format(nat,len(bulk)))
            csize=bulk.get_cell()                                                                                                         
            totalsol=Atoms(cell=csize, pbc=True)
            totalsol.extend(indi)
            totalsol.extend(bulk)
            for sym,c,m,u in Optimizer.atomlist:
                nc=len([atm for atm in totalsol if atm.symbol==sym])
                STR+='Defect configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n'
        elif Optimizer.structure=='Surface':
            totalsol=Atoms()
            totalsol.extend(indiv)
            nat=indiv.get_number_of_atoms()
            totalsol.extend(individ.bulki)
            if debug:
                logger.info('Extending surface structure to include bulk len(r1+r2)={0} len(bulk)={1}'.format(nat,len(individ.bulki)))
            for sym,c,m,u in Optimizer.atomlist:
                nc=len([atm for atm in totalsol if atm.symbol==sym])
                STR+='Surface-Bulk configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n'
            cell=numpy.maximum.reduce(indiv.get_cell())
            totalsol.set_cell([cell[0],cell[1],500])
            totalsol.set_pbc([True,True,False])
        elif Optimizer.structure=='Cluster':
           # logger.info('M:')
            totalsol = indiv.copy()
            nat = len(totalsol)
            if debug:
                logger.info('Extending cluster with {0} atoms to center of evaluation box of size {1}'.format(nat,Optimizer.large_box_size))
            origcell = indiv.get_cell()
            #print 'rank, eval_energy.cell',rank,origcell
            if Optimizer.forcing != 'RelaxBox':
               totalsol.set_cell([Optimizer.large_box_size,Optimizer.large_box_size,Optimizer.large_box_size])
               totalsol.translate([Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0])
           # logger.info('M: set cell')
        elif Optimizer.structure=='Crystal':
            totalsol = indiv.copy()
            nat = len(totalsol)
        else:
            print 'WARNING: In EvalEnergy. Optimizer.structure not recognized'
            logger.warning('Optimizer.structure not recognized')
        
        # Check for atoms that are too close or out of constrained location
        if Optimizer.constrain_position:
            if Optimizer.structure=='Defect':
                if debug:
                    logger.info('Constraining positions of defect')
                totalsol, stro = constrain_positions(totalsol, Optimizer.solidbulk, Optimizer.sf)
                if debug:
                    logger.info(stro)
                STR+=str0
        min_len=0.7
        if not Optimizer.fixed_region:
            if debug:
                logger.info('Running check minimum distance')
            # logger.info('M:check dist')
            totalsol, STR = check_min_dist(Optimizer, totalsol, Optimizer.structure, nat, min_len, STR)
            if debug:
                write_xyz(Optimizer.debugfile,totalsol,'After minlength check')
                Optimizer.debugfile.flush()
                logger.debug('Writing individual after checking minimum length')
        
        # Set calculator to use to get forces/energies
        if Optimizer.parallel:
           # logger.info('M:start calculator')
            calc = setup_calculator(Optimizer)
            if Optimizer.fixed_region:
                if debug:
                    logger.info('Setting up fixed region calculator')
                pms=copy.deepcopy(calc.parameters)
                try:
                    pms['mass'][len(pms['mass'])-1] += '\ngroup RO id >= {0}\nfix freeze RO setforce 0.0 0.0 0.0\n'.format(nat)
                except KeyError:
                    pms['pair_coeff'][0] += '\ngroup RO id >= {0}\nfix freeze RO setforce 0.0 0.0 0.0\n'.format(nat)
                calc = LAMMPS(parameters=pms, files=calc.files, keep_tmp_files=calc.keep_tmp_files, tmp_dir=calc.tmp_dir)
                lmin = copy.copy(Optimizer.lammps_min)
                if debug:
                    logger.info('Setting up no local minimization calculator')
                Optimizer.lammps_min = None
                Optimizer.static_calc = setup_calculator(Optimizer)
                Optimizer.lammps_min = lmin
        else:
            calc=Optimizer.calc
        totalsol.set_calculator(calc)
        totalsol.set_pbc(True)
        
        # Perform Energy Minimization
        if not Optimizer.parallel:
            if debug: 
                write_xyz(Optimizer.debugfile,totalsol,'Individual sent to Energy Minimizer')
                logger.debug('Writing structure sent to energy minimizer')
        try:
            cwd = os.getcwd()
            if Optimizer.ase_min == True:
                if debug:
                    logger.info('Running ASE minimizer')
                if Optimizer.calc_method=='LennardJones':
                    logger.warn('Must run ase LJ calculator with pbc=False')
                    totalsol.set_pbc(False)
                totalsol, energy, pressure, volume, STR = run_ase_min(totalsol, Optimizer.ase_min_fmax, Optimizer.ase_min_maxsteps, Optimizer.fitness_scheme, STR)
            else:
                if debug:
                    logger.info('Running local energy calculator')
                if Optimizer.fixed_region:
                    totalsol, pea, energy, pressure, volume, STR = run_energy_eval(totalsol, Optimizer.calc_method, Optimizer.fixed_region, Optimizer.fitness_scheme, STR, Optimizer.static_calc)
                else:
                  #  logger.info('M:start run_energy_eval')
                    totalsol, pea, energy, pressure, volume, STR = run_energy_eval(totalsol, Optimizer.calc_method, False, Optimizer.fitness_scheme, STR)
                    logger.info('M:finish run_energy_eval, energy = {0} @ rank ={1}'.format(energy,rank))
        except Exception, e:
            logger.critical('Error in energy evaluation: {0}'.format(e), exc_info=True)
            path = os.path.join(cwd,'TroubledLammps')
            if not os.path.exists(path):
                os.mkdir(path)
            #Copy files over
            shutil.copyfile(calc.trajfile,os.path.join(path,os.path.basename(calc.trajfile)))
            shutil.copyfile(calc.infile,os.path.join(path,os.path.basename(calc.infile)))
            shutil.copyfile(calc.logfile,os.path.join(path,os.path.basename(calc.logfile)))
            shutil.copyfile(calc.datafile,os.path.join(path,os.path.basename(calc.datafile)))
            raise RuntimeError('{0}:{1}'.format(Exception,e))
        if not Optimizer.parallel:
            if debug:
                write_xyz(Optimizer.debugfile,totalsol,'Individual after Energy Minimization')
                Optimizer.debugfile.flush()
                logger.debug('Writing structure recieved from energy minimizer')
       
        # Separate structures into distinct pieces
        if Optimizer.structure=='Defect':
            if Optimizer.fixed_region==True or Optimizer.finddefects==False:
                if debug:
                    logger.info('Identifying atoms in defect structure based on ID')
                individ[0]=totalsol[0:nat]
                bul=totalsol[(nat):len(totalsol)]
                individ[0].set_cell(csize)
            else:
                if debug:
                    logger.info('Applying find defects scheme to identify R1 and R2 for Defect')
                if 'FD' in Optimizer.debug:
                    outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=Optimizer.debugfile)
                else:
                    outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=False)
                individ[0]=outt[0]
                bul=outt[1]
                individ.vacancies = outt[2]
                individ.swaps = outt[3]
                STR += outt[4]
            indiv=individ[0]
        elif Optimizer.structure=='Surface':
            if debug:
                logger.info('Finding surface top layer')
            top,bul=find_top_layer(totalsol,Optimizer.surftopthick)
            indiv=top.copy()
            individ[0]=top.copy()
            bul = Atoms()
        elif Optimizer.structure=='Crystal':
            if debug:
                logger.info('Checking crystal cell type')
            celltype = check_cell_type(totalsol)
            STR+='Cell structure = {0}\n'.format(celltype)
            bul = Atoms()
            individ[0] = totalsol.copy()
        elif Optimizer.structure=='Cluster':
            volume = get_cluster_volume(totalsol)
            bul = Atoms()
            if debug:
                logger.info('Translating cluster back to smaller box size location')
            if Optimizer.forcing != 'RelaxBox':
               totalsol.translate([-Optimizer.large_box_size/2.0,-Optimizer.large_box_size/2.0,-Optimizer.large_box_size/2.0])
               totalsol.set_cell(origcell)
            individ[0] = totalsol.copy()
        
        # Add concentration energy dependence
        if Optimizer.forcing=='energy_bias':
            if debug:
                logger.info('Applying energy bias for atoms with different number of atoms of type than in atomlist')
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                n[i]=len([inds for inds in totalsol if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=abs(n[i]-Optimizer.atomlist[i][1])
            factor=sum(n)**3
            energy=(energy+factor)/totalsol.get_number_of_atoms()
            STR+='Energy with Bias = {0}\n'.format(energy)
        elif Optimizer.forcing=='chem_pot':
            if debug:
                logger.info('Applying chemical potential bias for atoms with different number of atoms of type than in atomlist')
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                n[i]=len([inds for inds in totalsol if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=n[i]*Optimizer.atomlist[i][3]
            factor=sum(n)
            energy=(energy+factor)/totalsol.get_number_of_atoms()
            STR+='Energy with Chemical Potential = {0}\n'.format(energy)

        individ.energy=energy
        individ.buli=bul
        individ.pressure=pressure
        individ.volume=volume

        #Add pealist to include atom index based on sorted PE. 
        logger.info('before sort{0}'.format(individ.energy))
        sort_pealist(Optimizer,individ,pea)
        energy = individ.energy
        logger.info('after sort {0}'.format(individ.energy))
        if Optimizer.fingerprinting:
            if debug:
                logger.info('Identifying fingerprint of new structure')
            individ.fingerprint=get_fingerprint(Optimizer,individ,Optimizer.fpbin,Optimizer.fpcutoff)
        if Optimizer.parallel:
            calc.clean()
            signal = 'Evaluated individual {0} on {1}\n'.format(individ.index,rank)
            signal +=STR
        else:
            signal=STR

    return energy, bul, individ, signal

def sort_pealist(Optimizer,individ,pea):
    logger = logging.getLogger(Optimizer.loggername)

    #Add pealist to include atom index based on sorted PE. 
    syms = [sym for sym,c,m,u in Optimizer.atomlist] 
    numatom = [c for sym,c,m,u in Optimizer.atomlist]
    peatom = [u for sym,c,m,u in Optimizer.atomlist]

    hpealist = [] 
    lpealist = []         
    for j in range(len(syms)) :
        sym = syms[j]
        pelist = []
        peindexlist = []
        for i in range(len(pea)) :
          if individ[0][i].symbol == sym:
            if pea[i][0] < peatom[j] - 5.0 or pea[i][0] > peatom[j] + 5.0 :
               individ.energy = 10000
               message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(individ.history_index)
               logger.warn(message)
            pelist.append(pea[i][0])
            peindexlist.append(i)     
        pearray = numpy.asarray(pelist) 
        pearray_sorted = pearray.argsort()     
        for i in range(1,11):
           hpealist.append(peindexlist[pearray_sorted[-i]])
        for i in range(0,10):
           lpealist.append(peindexlist[pearray_sorted[i]])
            
        individ.hpealist = hpealist
        individ.lpealist = lpealist
    return 


def constrain_positions(indiv, bulk, sf):
    STR=''
    ts = indiv.copy()
    indc,indb,vacant,swap,stro = find_defects(ts,bulk,0)
    sbulk = bulk.copy()
    bcom = sbulk.get_center_of_mass()
    #totalsol.translate(-bulkcom)
    #indc.translate(-bulkcom)
    #totalsol.append(Atom(position=[0,0,0]))
    #             for one in indc:
    #                 index = [atm.index for atm in totalsol if atm.position[0]==one.position[0] and atm.position[1]==one.position[1] and atm.position[2]==one.position[2]][0]
    #                 if totalsol.get_distance(-1,index) > Optimizer.sf:
    #                     r = random.random()
    #                     totalsol.set_distance(-1,index,Optimizer.sf*r,fix=0)
    #             totalsol.pop()
    #             totalsol.translate(bulkcom)
    com = indc.get_center_of_mass()
    dist = (sum((bcom[i] - com[i])**2 for i in range(3)))**0.5
    if dist > sf:
        STR+='Shifting structure to within region\n'
        r = random.random()*sf
        comv = numpy.linalg.norm(com)
        ncom = [one*r/comv for one in com]
        trans = [ncom[i]-com[i] for i in range(3)]
        indices = []
        for one in indc:
            id = [atm.index for atm in totalsol if atm.position[0]==one.position[0] and atm.position[1]==one.position[1] and atm.position[2]==one.position[2]][0]
            ts[id].position += trans
    return ts, STR

def check_min_dist(Optimizer, totalsol, type='Defect', nat=None, min_len=0.7, STR=''):
    if type=='Defect' or type=='Crystal' or type=='Surface':
        if nat==None:
            nat=len(totalsol)
        cutoffs=[2.0 for one in totalsol]
        nl=NeighborList(cutoffs,bothways=True,self_interaction=False)
        nl.update(totalsol)
        for one in totalsol[0:nat]:
            nbatoms=Atoms()
            nbatoms.append(one)
            indices, offsets=nl.get_neighbors(one.index)
            for index, d in zip(indices,offsets):
                index = int(index)
                sym=totalsol[index].symbol
                pos=totalsol[index].position + numpy.dot(d,totalsol.get_cell())
                at=Atom(symbol=sym,position=pos)
                nbatoms.append(at)
            while True:
                dflag=False
                for i in range(1,len(nbatoms)):
                    d=nbatoms.get_distance(0,i)
                    if d < min_len:
                        nbatoms.set_distance(0,i,min_len+.01,fix=0.5)
                        STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
                        dflag=True
                if dflag==False:
                    break
            for i in range(len(indices)):
                totalsol[indices[i]].position=nbatoms[i+1].position
            totalsol[one.index].position=nbatoms[0].position
            nl.update(totalsol)
    elif type=='Cluster':
        rank = MPI.COMM_WORLD.Get_rank()
        logger = logging.getLogger(Optimizer.loggername)
        R = totalsol.arrays['positions']
        tol = 0.01
        epsilon = 0.05
        fix = 0.5
        if Optimizer.forcing == 'EllipoidShape' or Optimizer.forcing == 'FreeNatom': 
          com = totalsol.get_center_of_mass()       
          cmax = numpy.maximum.reduce(R)
          cmin = numpy.minimum.reduce(R)
          rmax= (cmax-cmin)/2.0 
          if Optimizer.forcing == 'FreeNatom':
             rcutoff = 44.0
             cutoff = [44.0,44.0,20.0]        
          else:
             rcutoff = 11.0
             cutoff = [12.0,12.0,12.0]        
          rcutoff = 44.0
          cutoff = [44.0,44.0,20.0]        
          #check if atoms are isolated outside of cluster
          cutoffs=[3.0 for one in totalsol]
          nl=NeighborList(cutoffs,bothways=True,self_interaction=False)
          nl.update(totalsol)
          for i in range(len(totalsol)):
             indices, offsets=nl.get_neighbors(i)
             D = R[i]-com
             radius = (numpy.dot(D,D))**0.5 #numpy.linalg.norm(D)
             if len(indices) < 12 or radius > rcutoff :
                # logger.info('M:move atoms back when indice {0} or radius {1}'.format(len(indices),radius))
                # R[i] = [com[j] + D[j]/radius*rcutoff for j in range(3)]
                theta=math.radians(random.uniform(0,360))
                phi=math.radians(random.uniform(0,180))
                R[i][0] = com[0] + (rmax[0]+2.5)*math.sin(theta)*math.cos(phi) #allow atoms expend by 2.5 ang
                R[i][1] = com[1] + (rmax[1]+2.5)*math.sin(theta)*math.sin(phi)
                R[i][2] = com[2] + rmax[2]*math.cos(theta)
                # logger.info('M:move atoms back new pos {0} {1} {2}'.format(rmax[0]*math.sin(theta)*math.cos(phi),rmax[1]*math.sin(theta)*math.sin(phi),rmax[2]*math.cos(theta)))
               
          # check if atoms are within cluster region 
          for i in range(0,len(totalsol)):
            # D = R[i]-com
            for j in range(3):                 
               if D[j] > cutoff[j] or D[j] < -cutoff[j]:
         #         logger.info('M:before move R {0} & com {1}'.format(R[i][j],com[j]))
                  #if rank==0:
                  #   print "before:",j,R[i][j],com[j] 
                  R[i][j] = com[j]+numpy.sign(D[j])*cutoff[j]*random.random()
         #         logger.info('M:after move R {0} '.format(R[i][j]))
                  #if rank==0:
                  #   print "after:",R[i][j]
              #    STR+='--- WARNING: Atoms too far along x-y (>44A) - Implement Move ---\n'          
                  D = R[i]-com
        #    radius = (numpy.dot(D,D))**0.5 #numpy.linalg.norm(D)
             #  STR+='--- WARNING: Atoms too far (>56A) - Implement Move ---\n'          
              

        closelist = numpy.arange(len(totalsol))
        iter = 0
        while len(closelist) > 0 and iter<2:
          iter+=1 
         # checklist = numpy.copy(closelist)
          closelist = []  
          dist=scipy.spatial.distance.cdist(R,R)       
          numpy.fill_diagonal(dist,1.0)
          smalldist = numpy.where(dist < min_len-tol)
         # for i in checklist:
            # for j in range(i+1,len(totalsol)):
         #    if len(checklist) == len(totalsol):
         #       jstart = i+1
         #    else:
         #       jstart = 0
         #    for j in range(jstart,len(totalsol)):
         #       if i != j and dist[i][j] < min_len:
             #       d=totalsol.get_distance(i,j)
             #       if d < min_len:
             #           totalsol.set_distance(i,j,min_len,fix=0.5)
                    # d = (D[0]*D[0]+D[1]*D[1]+D[2]*D[2])**0.5
          for ind in range(len(smalldist[0])):
                   i = smalldist[0][ind]
                   j = smalldist[1][ind]
                   if i < j and dist[i][j] < min_len-tol:   
                        closelist.append(i)
                        closelist.append(j)
                        if dist[i][j] > epsilon:
                      	  x = 1.0 - min_len / dist[i][j]
                          D = R[j]-R[i]
                         # print "node:",rank,"x",x,R[i],R[j],D, dist[i][j]
                       	  R[i] += (x * fix) * D
                          R[j] -= (x * (1.0 - fix)) * D
                        else:
                          R[i] += [0.2, 0.0, 0.0]
                          R[j] -= [0.2, 0.0, 0.0] 
                        R2P = [R[i],R[j]]
                        dist2P=scipy.spatial.distance.cdist(R2P,R)       
                        dist[i] = dist2P[0]
                        dist[j] = dist2P[1]
                        for k in range(len(R)):
                            dist[k][i] = dist[i][k]
                            dist[k][j] = dist[j][k]
                      #  STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
          closelist=list(set(closelist))
          closelist.sort()
          if len(closelist) != 0: 
             logger.info('M:iter {0}, closelist size {1}'.format(iter,len(closelist)))
         #    print "rank", rank, closelist
    else:
        print 'WARNING: In Check_Min_Dist in EvalEnergy: Structure Type not recognized'
    return totalsol, STR

def get_cluster_volume(cluster):
    max_pos = numpy.maximum.reduce(cluster.get_positions())
    min_pos = numpy.minimum.reduce(cluster.get_positions())
    diff_pos = [max_pos[i]-min_pos[i] for i in range(3)]
    vol = diff_pos[0]*diff_pos[1]*diff_pos[2]
    return vol

def run_ase_min(totalsol, fmax=0.01, mxstep=1000, fitscheme='totalenfit', STR=''):
    try:
        dyn=BFGS(totalsol)
        dyn.run(fmax=fmax, steps=mxstep)
    except OverflowError:
        STR+='--- Error: Infinite Energy Calculated - Implement Random shake---\n'
        totalsol.rattle(stdev=0.3)
        dyn=BFGS(totalsol)
        dyn.run(fmax=fmax, steps=mxstep)
    except numpy.linalg.linalg.LinAlgError:
        STR+='--- Error: Singular Matrix - Implement Random shake ---\n'
        totalsol.rattle(stdev=0.2)
        dyn=BFGS(totalsol)
        dyn.run(fmax=fmax, steps=mxstep)
    # Get Energy of Minimized Structure
    en=totalsol.get_potential_energy()
    if fitscheme == 'enthalpyfit':
        pressure=totalsol.get_isotropic_pressure(totalsol.get_stress())
    else:
        pressure=0
    volume = totalsol.get_volume()
    energy=en
    return totalsol, energy, pressure, volume, STR

def run_energy_eval(totalsol, calc_method='LAMMPS', fx_region=False, fit_scheme='totalenfit', STR='', static_calc=None):
    if calc_method=='VASP':
        en=totalsol.get_potential_energy()
        calcb=Vasp(restart=True)
        totalsol=calcb.get_atoms()
        stress=calcb.read_stress()
    else:
        totcop = totalsol.copy()
        OUT = totalsol.calc.calculate(totalsol)
        totalsol = OUT['atoms']
        pea = OUT['pea']
       # M: test
       # if MPI.COMM_WORLD.Get_rank() == 0:
       #   for i in range(len(totalsol)) :
       #      print i,totalsol[i].symbol,totalsol[i].position, pea[i]
        totalsol.set_pbc(True)
        if fx_region:
            STR+='Energy of fixed region calc = {0}\n'.format(OUT['thermo'][-1]['pe'])
            totalsol.set_calculator(static_calc)
            OUT=totalsol.calc.calculate(totalsol)
            totalsol=OUT['atoms']
            totalsol.set_pbc(True)
            STR+='Energy of static calc = {0}\n'.format(OUT['thermo'][-1]['pe'])
        en=OUT['thermo'][-1]['pe']
        stress=numpy.array([OUT['thermo'][-1][i] for i in ('pxx','pyy','pzz','pyz','pxz','pxy')])*(-1e-4*GPa)
    if fit_scheme == 'enthalpyfit':
        pressure = totalsol.get_isotropic_pressure(stress)
    else:
        pressure = 0
    volume = totalsol.get_volume()
    energy=en
    STR+='Energy per atom = {0}\n'.format(energy/len(totalsol))
    return totalsol, pea, energy, pressure, volume, STR
