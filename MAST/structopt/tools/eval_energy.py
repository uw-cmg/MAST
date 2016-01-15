import os
from ase import Atom, Atoms
from MAST.structopt.inp_out.write_xyz import write_xyz
from MAST.structopt.tools.setup_calculator import setup_calculator
from MAST.structopt.tools.find_defects import find_defects
from MAST.structopt.tools.check_cell_type import check_cell_type
from MAST.structopt.fingerprinting import get_fingerprint
from MAST.structopt.tools.lammps import LAMMPS
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.units import GPa
from ase.calculators.neighborlist import NeighborList
import numpy
try:
    from mpi4py import MPI
except ImportError:
    pass
import logging
import pdb
import shutil

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
            totalsol = indiv.copy()
            nat = len(totalsol)
            if debug:
                logger.info('Extending cluster with {0} atoms to center of evaluation box of size {1}'.format(nat,Optimizer.large_box_size))
            origcell = indiv.get_cell()
            totalsol.set_cell([Optimizer.large_box_size,Optimizer.large_box_size,Optimizer.large_box_size])
            totalsol.translate([Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0])
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
            totalsol, STR = check_min_dist(totalsol, Optimizer.structure, nat, min_len, STR)
            if debug:
                write_xyz(Optimizer.debugfile,totalsol,'After minlength check')
                Optimizer.debugfile.flush()
                logger.debug('Writing individual after checking minimum length')
        
        # Set calculator to use to get forces/energies
        if Optimizer.parallel:
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
                    totalsol, energy, pressure, volume, STR = run_energy_eval(totalsol, Optimizer.calc_method, Optimizer.fixed_region, Optimizer.fitness_scheme, STR, Optimizer.static_calc)
                else:
                    totalsol, energy, pressure, volume, STR = run_energy_eval(totalsol, Optimizer.calc_method, False, Optimizer.fitness_scheme, STR)
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

def check_min_dist(totalsol, type='Defect', nat=None, min_len=0.7, STR=''):
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
        for i in range(len(totalsol)):
            for j in range(len(totalsol)):
                if i != j:
                    d=totalsol.get_distance(i,j)
                    if d < min_len:
                        totalsol.set_distance(i,j,min_len,fix=0.5)
                        STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
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
    return totalsol, energy, pressure, volume, STR
