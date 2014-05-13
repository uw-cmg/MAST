import os
from ase import Atom, Atoms
from MAST.structopt.io.write_xyz import write_xyz
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
import pdb

def eval_energy(input):
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
    if input[0]==None:
        energy=0
        bul=0
        individ=0
        rank = MPI.COMM_WORLD.Get_rank()
        signal='Evaluated none individual on '+repr(rank)+'\n'
    else:
        [Optimizer, individ]=input
    if Optimizer.calc_method=='MAST':
        energy = individ.energy
        bul = individ.energy
        signal = 'Recieved MAST structure\n'
    else:
        if Optimizer.parallel: rank = MPI.COMM_WORLD.Get_rank()
        if not Optimizer.genealogy:
            STR='----Individual ' + str(individ.index)+ ' Optimization----\n'
        else:
            STR='----Individual ' + str(individ.history_index)+ ' Optimization----\n'
        indiv=individ[0]
        if 'EE' in Optimizer.debug:
            debug = True
        else:
            debug = False
        if debug: 
            write_xyz(Optimizer.debugfile,indiv,'Recieved by eval_energy')
            Optimizer.debugfile.flush()
        if Optimizer.structure=='Defect':
            indi=indiv.copy()
            if Optimizer.alloy==True:
                bulk=individ.bulki
            else:
                bulk=individ.bulko
            nat=indi.get_number_of_atoms()
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
            for sym,c,m,u in Optimizer.atomlist:
                nc=len([atm for atm in totalsol if atm.symbol==sym])
                STR+='Surface-Bulk configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n'
            cell=numpy.maximum.reduce(indiv.get_cell())
            totalsol.set_cell([cell[0],cell[1],500])
            totalsol.set_pbc([True,True,False])
    
        if Optimizer.constrain_position:
            ts = totalsol.copy()
            indc,indb,vacant,swap,stro = find_defects(ts,Optimizer.solidbulk,0)
            sbulk = Optimizer.solidbulk.copy()
            bcom = sbulk.get_center_of_mass()
            #totalsol.translate(-bulkcom)
            #indc.translate(-bulkcom)
            #totalsol.append(Atom(position=[0,0,0]))
    # 			for one in indc:
    # 				index = [atm.index for atm in totalsol if atm.position[0]==one.position[0] and atm.position[1]==one.position[1] and atm.position[2]==one.position[2]][0]
    # 				if totalsol.get_distance(-1,index) > Optimizer.sf:
    # 					r = random.random()
    # 					totalsol.set_distance(-1,index,Optimizer.sf*r,fix=0)
    # 			totalsol.pop()
    # 			totalsol.translate(bulkcom)
            com = indc.get_center_of_mass()
            dist = (sum((bcom[i] - com[i])**2 for i in range(3)))**0.5
            if dist > Optimizer.sf:
                STR+='Shifting structure to within region\n'
                r = random.random()*Optimizer.sf
                comv = numpy.linalg.norm(com)
                ncom = [one*r/comv for one in com]
                trans = [ncom[i]-com[i] for i in range(3)]
                indices = []
                for one in indc:
                    id = [atm.index for atm in totalsol if atm.position[0]==one.position[0] and atm.position[1]==one.position[1] and atm.position[2]==one.position[2]][0]
                    totalsol[id].position += trans
    
        # Check for atoms that are too close
        min_len=0.7
        #pdb.set_trace()
        if not Optimizer.fixed_region:
            if Optimizer.structure=='Defect' or Optimizer.structure=='Surface':
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
                if debug:
                    write_xyz(Optimizer.debugfile,totalsol,'After minlength check')
                    Optimizer.debugfile.flush()
            else:
                for i in range(len(indiv)):
                    for j in range(len(indiv)):
                        if i != j:
                            d=indiv.get_distance(i,j)
                            if d < min_len:
                                indiv.set_distance(i,j,min_len,fix=0.5)
                                STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
                if debug:
                    write_xyz(Optimizer.debugfile,indiv,'After minlength check')
                    Optimizer.debugfile.flush()
    
        # Set calculator to use to get forces/energies
        if Optimizer.parallel:
            calc = setup_calculator(Optimizer)
            if Optimizer.fixed_region:
                pms=copy.deepcopy(calc.parameters)
                try:
                    pms['mass'][len(pms['mass'])-1] += '\ngroup RO id >= '+repr(nat)+'\nfix freeze RO setforce 0.0 0.0 0.0\n'
                except KeyError:
                    pms['pair_coeff'][0] += '\ngroup RO id >= '+repr(nat)+'\nfix freeze RO setforce 0.0 0.0 0.0\n'
                calc = LAMMPS(parameters=pms, files=calc.files, keep_tmp_files=calc.keep_tmp_files, tmp_dir=calc.tmp_dir)
                lmin = copy.copy(Optimizer.lammps_min)
                Optimizer.lammps_min = None
                Optimizer.static_calc = setup_calculator(Optimizer)
                Optimizer.lammps_min = lmin
        else:
            calc=Optimizer.calc
        if Optimizer.structure=='Defect' or Optimizer.structure=='Surface':
            totalsol.set_calculator(calc)
            totalsol.set_pbc(True)
        else:
            indiv.set_calculator(calc)
            indiv.set_pbc(True)	#Current bug in ASE optimizer-Lammps prevents pbc=false 
            if Optimizer.structure=='Cluster':
                indiv.set_cell([500,500,500])
                indiv.translate([250,250,250])
    
        cwd=os.getcwd()
        # Perform Energy Minimization
        if not Optimizer.parallel:
            Optimizer.output.flush()
        if Optimizer.ase_min == True:
            try:
                if Optimizer.structure=='Defect' or Optimizer.structure=='Surface':
                    dyn=BFGS(totalsol)
                else:
                    dyn=BFGS(indiv)
                dyn.run(fmax=Optimizer.ase_min_fmax, steps=Optimizer.ase_min_maxsteps)
            except OverflowError:
                STR+='--- Error: Infinite Energy Calculated - Implement Random ---\n'
                box=Atoms()
                indiv=gen_pop_box(Optimizer.natoms, Optimizer.atomlist, Optimizer.size)
                indiv.set_calculator(calc)
                dyn=BFGS(indiv)
                dyn.run(fmax=fmax, steps=steps)
            except numpy.linalg.linalg.LinAlgError:
                STR+='--- Error: Singular Matrix - Implement Random ---\n'
                indiv=gen_pop_box(Optimizer.natoms, Optimizer.atomlist, Optimizer.size)
                indiv.set_calculator(calc)
                dyn=BFGS(indiv)
                dyn.run(fmax=fmax, steps=steps)
            # Get Energy of Minimized Structure
            if Optimizer.structure=='Defect' or Optimizer.structure=='Surface':
                en=totalsol.get_potential_energy()
                #force=numpy.maximum.reduce(abs(totalsol.get_forces()))
                if Optimizer.fitness_scheme == 'enthalpyfit':
                    pressure=totalsol.get_isotropic_pressure(totalsol.get_stress())
                    cell_max=numpy.maximum.reduce(totalsol.get_positions())
                    cell_min=numpy.minimum.reduce(totalsol.get_positions())
                    cell=cell_max-cell_min
                    volume=cell[0]*cell[1]*cell[2]
                else:
                    pressure=0
                    volume=0
                na=totalsol.get_number_of_atoms()
                ena=en/na
                energy=en
                individ[0]=totalsol[0:nat]
                bul=totalsol[(nat):len(totalsol)]
                STR+='Number of positions = '+repr(len(bul)+len(individ[0]))+'\n'
                individ[0].set_cell(csize)
                indiv=individ[0]
            else:
                en=indiv.get_potential_energy()
                if Optimizer.fitness_scheme == 'enthalpyfit':
                    pressure=indiv.get_isotropic_pressure(indiv.get_stress())
                    cell_max=numpy.maximum.reduce(indiv.get_positions())
                    cell_min=numpy.minimum.reduce(indiv.get_positions())
                    cell=cell_max-cell_min
                    volume=cell[0]*cell[1]*cell[2]
                else: 
                    pressure=0
                    volume=0
                na=indiv.get_number_of_atoms()
                ena=en/na
                energy=ena
                individ[0]=indiv
                bul=0
        else:
            if Optimizer.structure=='Defect' or Optimizer.structure=='Surface':
                if Optimizer.calc_method=='VASP':
                    en=totalsol.get_potential_energy()
                    calcb=Vasp(restart=True)
                    totalsol=calcb.get_atoms()
                    stress=calcb.read_stress()
                else:
                    try:
                        totcop=totalsol.copy()
                        if debug: write_xyz(Optimizer.debugfile,totcop,'Individual sent to lammps')
                        OUT=totalsol.calc.calculate(totalsol)
                        totalsol=OUT['atoms']
                        totalsol.set_pbc(True)
                        if Optimizer.fixed_region:
                            if debug:
                                print 'Energy of fixed region calc = ', OUT['thermo'][-1]['pe']
                            totalsol.set_calculator(Optimizer.static_calc)
                            OUT=totalsol.calc.calculate(totalsol)
                            totalsol=OUT['atoms']
                            totalsol.set_pbc(True)
                            if debug:
                                print 'Energy of static calc = ', OUT['thermo'][-1]['pe']
                        en=OUT['thermo'][-1]['pe']
                        stress=numpy.array([OUT['thermo'][-1][i] for i in ('pxx','pyy','pzz','pyz','pxz','pxy')])*(-1e-4*GPa)
                        #force=numpy.maximum.reduce(abs(totalsol.get_forces()))
                        if debug:
                            write_xyz(Optimizer.debugfile,totalsol,'After Lammps Minimization')
                            Optimizer.debugfile.flush()
                    except Exception, e:
                        os.chdir(cwd)
                        STR+='WARNING: Exception during energy eval:\n'+repr(e)+'\n'
                        f=open('problem-structures.xyz','a')
                        write_xyz(f,totcop,data='Starting structure hindex='+individ.history_index)
                        write_xyz(f,totalsol,data='Lammps Min structure')
                        en=10
                        stress=0
                        f.close()
                if Optimizer.fitness_scheme == 'enthalpyfit':
                    pressure=totalsol.get_isotropic_pressure(stress)
                    cell_max=numpy.maximum.reduce(totalsol.get_positions())
                    cell_min=numpy.minimum.reduce(totalsol.get_positions())
                    cell=cell_max-cell_min
                    volume=cell[0]*cell[1]*cell[2]
                else:
                    pressure=totalsol.get_isotropic_pressure(stress)
                    volume=0
                na=totalsol.get_number_of_atoms()
                ena=en/na
                energy=en
                if Optimizer.structure=='Defect':
                    if Optimizer.fixed_region==True or Optimizer.finddefects==False:
                        individ[0]=totalsol[0:nat]
                        bul=totalsol[(nat):len(totalsol)]
                        individ[0].set_cell(csize)
                    else:
                        if 'FI' in Optimizer.debug:
                            outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=Optimizer.debugfile)
                        else:
                            outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=False)
                        individ[0]=outt[0]
                        bul=outt[1]
                        individ.vacancies = outt[2]
                        individ.swaps = outt[3]
                        STR += outt[4]
                    indiv=individ[0]
                else:
                    top,bul=find_top_layer(totalsol,Optimizer.surftopthick)
                    indiv=top.copy()
                    individ[0]=top.copy()
            else:
                if Optimizer.calc_method=='VASP':
                    en=totalsol.get_potential_energy()
                    calcb=Vasp(restart=True)
                    totalsol=calcb.get_atoms()
                    stress=calcb.read_stress()
                else:
                    OUT=indiv.calc.calculate(indiv)
                    en=OUT['thermo'][-1]['pe']
                    #indiv.set_positions(OUT['atoms'].get_positions())
                    #indiv.set_cell(OUT['atoms'].get_cell())
                    indiv=OUT['atoms']
                    indiv.set_pbc(True)
                    stress=numpy.array([OUT['thermo'][-1][i] for i in ('pxx','pyy','pzz','pyz','pxz','pxy')])*(-1e-4*GPa)
                if Optimizer.fitness_scheme == 'enthalpyfit':
                    pressure=indiv.get_isotropic_pressure(stress)
                    cell_max=numpy.maximum.reduce(indiv.get_positions())
                    cell_min=numpy.minimum.reduce(indiv.get_positions())
                    cell=cell_max-cell_min
                    volume=cell[0]*cell[1]*cell[2]
                else: 
                    pressure=indiv.get_isotropic_pressure(stress)
                    volume=0
                na=indiv.get_number_of_atoms()
                ena=en/na
                energy=en
                individ[0]=indiv
                bul=0
        STR+='EnergypAtm = '+repr(ena)+'\n'
        if Optimizer.structure=='Crystal':
            STR+='Cell structure = '+repr(check_cell_type(indiv))+'\n'
    
        # Add concentration energy dependence
        if Optimizer.forcing=='Energy_bias':
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                n[i]=len([inds for inds in indiv if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=abs(n[i]-Optimizer.atomlist[i][1])
            factor=sum(n)**3
            energy=(en+factor)/na
        elif Optimizer.forcing=='Chem_pot':
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                if Optimizer.structure=='Defect':
                    n[i]=len([inds for inds in totalsol if inds.symbol==Optimizer.atomlist[i][0]])
                else:
                    n[i]=len([inds for inds in indiv if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=n[i]*Optimizer.atomlist[i][3]
            factor=sum(n)
            energy=(en+factor)/na
            STR+='Energy with Chemical Potential = '+repr(energy[0])+'\n'

        # Add explosion prevention protection
    # 	if 'prevent_explosions' in globals():
    # 		dist=[10]*len(indiv)
    # 		for i in range(len(indiv)):
    # 			for j in range(len(indiv)):
    # 				if i != j:
    # 					dist[j]=indiv.get_distance(i,j,mic=True)
    # 			if min(dist) > 3.5:
    # 				energy+=10
    
        #if Optimizer.structure=='Defect':
        #	individ.force=force
        individ.energy=energy
        individ.buli=bul
        individ.pressure=pressure
        individ.volume=volume
    
        if Optimizer.structure=='Cluster':
            indiv.translate([-250,-250,-250])
        if Optimizer.fingerprinting:
            individ.fingerprint=get_fingerprint(Optimizer,individ,Optimizer.fpbin,Optimizer.fpcutoff)
        if Optimizer.parallel:
            calc.clean()
            signal = 'Evaluated individual '+repr(individ.index)+' on '+repr(rank)+'\n'
            signal +=STR
        else:
            signal=STR

    return energy, bul, individ, signal

