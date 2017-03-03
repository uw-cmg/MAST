import sys
import os
import time
import random
import math
import pdb
import logging
from MAST.structopt import inp_out
from MAST.structopt import tools
from MAST.structopt import generate
from MAST.structopt import switches
from MAST.structopt.switches.predator_switch import predator_switch
from MAST.structopt import fingerprinting
from MAST.structopt import post_processing as pp
from ase import Atom, Atoms
try:
    from mpi4py import MPI
except ImportError:
    pass

class Optimizer():
    __version__  = 'StructOpt_v0.4'
    logger = None
    
    def __init__(self, input, uselogger=True):
        raise DeprecationWarning("StructOpt packaged with MAST is now deprecated. Please see https://github.com/uw-cmg/StructOpt_modular instead.")
        if input:
            parameters = inp_out.read_parameter_input(input, uselogger)
        else:
            parameters = inp_out.read_parameter_input({'atomlist':[('Xx',1,0,0)],'structure':'Cluster'}, uselogger)
        self.__dict__.update(parameters)
        if self.loggername:
            global logger
            logger = logging.getLogger(self.loggername)
        if self.restart_optimizer:
            try:
                rank = MPI.COMM_WORLD.Get_rank()
            except:
                rank = 0
            if rank==0:
                logger.info('restarting output')
                outdict = inp_out.restart_output(self)
                self.__dict__.update(outdict)
                logger.info('Loading individual files')
                poplist = []
                for indfile in self.population:
                    ind = inp_out.read_individual(indfile)
                    poplist.append(ind)
                self.population = poplist
                logger.info('Loading bests')
                bestlist = []
                for bestfile in self.BESTS:
                    ind = inp_out.read_individual(bestfile)
                    bestlist.append(ind)
                self.BESTS = bestlist
                self.restart = True
                if self.structure == 'Defect':
                    bulk = inp_out.read_xyz(self.solidfile)
                    bulk.set_pbc(True)
                    bulk.set_cell(self.solidcell)
                    self.solidbulk = bulk.copy()
                else:
                    self.solidbulk = None
        else:
            self.convergence = False
            self.generation = 0
            self.Runtimes = [time.time()]
            self.Evaluations = list()
            self.CXs = list()
            self.Muts = list()
            self.cxattempts = 0
            self.mutattempts = list()
            self.BESTS = list()
            self.genrep = 0
            self.minfit = 0
            self.convergence = False
            self.overrideconvergence = False
            self.population = list()
            self.calc = None
            self.static_calc = None
    
    def algorithm_initialize(self):
        global logger
        if self.restart_optimizer:
            logger.info('Successfully loaded optimizer from file')
        else:
            #Setup the output files
            logger.info('Initializing output for algorithm')
            outdict = inp_out.setup_output(self.filename, self.restart, self.nindiv,
                self.indiv_defect_write, self.genealogy, self.allenergyfile,
                self.fingerprinting, self.debug)
            self.__dict__.update(outdict)
            #Set starting convergence and generation
            logger.info('Initializing convergence, generation, and output monitoring stats')
            self.convergence = False
            self.overrideconvergence = False
            self.generation = 0
            #Prep output monitoring
            self.Runtimes = [time.time()]
            self.Evaluations = list()
            self.CXs = list()
            self.Muts = list()
            self.cxattempts = 0
            self.mutattempts = list()
            #Initialize random number generator
            random.seed(self.seed)
            #Initialize swap list variable
            if self.swaplist == None:
                #Use if concentrations in swaplist is unknown
                self.swaplist = [ [sym,c] for sym,c,m,u in self.atomlist ]
                swaplist = self.swaplist
            elif self.swaplist:
                swaplist = self.swaplist
            #Write the input parameters to the output file
            logger.info('Writing the input parameters to output file')
            inp_out.write_parameters(self)
    
    def algorithm_island(self):
        """Subprogram to run the GA in a Island Model Parallel Scheme"""
        global logger
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        if 'MA' in self.debug:
            debug = True
        else:
            debug = False
        if rank==0:
            #Initialize Random - Ensure each population has different random seed
            seeds=[]
            while len(seeds) < comm.Get_size():
                a=random.randint(0,100)
                if a not in seeds:
                    seeds.append(a)
        else:
            seeds=None
        self = comm.bcast(self,root=0)
        seed = comm.scatter(seeds,root=0)
        self.seed=seed
        random.seed(self.seed)
        self.algorithm_initialize()
        logger.info('Beginning main algorithm loop')
        logger.info('Random seed = {0}'.format(seed))
        MIGEVENT=0
        while True:
            pop = self.population
            #Migration
            if self.generation != 0:
                if float(self.generation)/float(self.migration_intervals)%1==0:
                    orcs=comm.allgather(self.overrideconvergence)
                    if True in orcs:
                        break
                    logger.info('Migrating structures')
                    self.output.write('-----Migration in generation '+repr(self.generation)+'-----\n')
                    nmig=int(self.nindiv*self.migration_percent)
                    migrantsind=[]
                    for i in range(nmig):
                        migrantsind.append(pop[i].duplicate())
                    self.output.write('Migrating '+repr(len(migrantsind))+' individuals\n')
                    migrantsind=[migrantsind]
                    migrants = comm.gather(migrantsind,root=0)
                    if rank==0:
                        self.output.write('Master recieved '+repr(len(migrants))+' sets of migrants\n')
                        MIGEVENT+=1
                        if MIGEVENT>len(migrants):
                            MIGEVENT-=len(migrants)
                        nmigrants=migrants[MIGEVENT:]+migrants[:MIGEVENT]
                    else:
                        nmigrants=None
                    nmigrants = comm.scatter(nmigrants,root=0)
                    for one in nmigrants[0]:
                        pop.append(one)
            offspring = self.generation_set(pop)
            # Identify the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if ind.fitness==0]
            #Evaluate the individuals with invalid fitness
            self.output.write('\n--Evaluate Structures--\n')
            poorlist=[]
            logger.info('Evaluating fitness of structures')
            for i in range(len(invalid_ind)):
                ind = invalid_ind[i]
                outs = switches.fitness_switch([self,ind])
                self.output.write(outs[1])
                invalid_ind[i] = outs[0]
            pop.extend(invalid_ind)
            pop = self.generation_eval(pop)
            if self.convergence:
                self.overrideconvergence=True
            self.write()
        
        logger.info('Run algorithm stats')
        end_signal = self.algorithm_stats(self.population)
        return end_signal    
    
    def algorithm_parallel1(self):
        """Subprogram for running parallel version of GA
        Requires MPI4PY"""
        global logger
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        if 'MA' in self.debug:
            debug = True
        else:
            debug = False
        if rank==0:
            self.algorithm_initialize()
            logger.info('Beginning main algorithm loop')
        #Begin main algorithm loop
        self.convergence = False
        convergence=False
        while not convergence:
            if rank==0:
                pop = self.population
                offspring = self.generation_set(pop)
                # Identify the individuals with an invalid fitness
                invalid_ind = [ind for ind in offspring if ind.fitness==0]
                #Evaluate the individuals with invalid fitness
                self.output.write('\n--Evaluate Structures--\n')
                ntimes=int(math.ceil(float(len(invalid_ind))/float(comm.Get_size())))
                nadd=int(ntimes*comm.Get_size()-len(invalid_ind))
                maplist=[[] for n in range(ntimes)]
                strt=0
                for i in range(len(maplist)):
                    maplist[i]=[[self,indi] for indi in invalid_ind[strt:comm.Get_size()+strt]]
                    strt+=comm.Get_size()
                for i in range(nadd):
                    maplist[len(maplist)-1].append([None,None])
#                 for i in range(ntimes):
#                     for j in range(
            else:
                ntimes=None
#             worker = MPI.COMM_SELF.Spawn(cmd, None, 5)
# 
#             n  = array('i', [100])
#             worker.Bcast([n,MPI.INT], root=MPI.ROOT)
# 
#             pi = array('d', [0.0])
#             worker.Reduce(sendbuf=None,
#                           recvbuf=[pi, MPI.DOUBLE],
#                           op=MPI.SUM, root=MPI.ROOT)
#             pi = pi[0]
# 
#             worker.Disconnect()
#             ntimes = comm.bcast(ntimes,root=0)
            outs=[]
            for i in range(ntimes):
                if rank==0:
                    one=maplist[i]
                else:
                    one=None
                ind =comm.scatter(one,root=0)
                out = switches.fitness_switch(ind)
                outt = comm.gather(out,root=0)
                if rank==0:
                    outs.extend(outt)
            if rank==0:
                for i in range(len(invalid_ind)):
                    invalid_ind[i] = outs[i][0]
                for i in range(len(outs)):
                    self.output.write(outs[i][1])
                pop.extend(invalid_ind)
                pop = self.generation_eval(pop)
                self.write()
            convergence =comm.bcast(self.convergence, root=0)
        
        if rank==0:
            logger.info('Run algorithm stats')
            end_signal = self.algorithm_stats(self.population)
        else:
            end_signal = None
        end_signal = comm.bcast(end_signal, root=0)
        return end_signal    
    
    def algorithm_parallel(self):
        """Subprogram for running parallel version of GA
        Requires MPI4PY"""
        global logger
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        if 'MA' in self.debug:
            debug = True
        else:
            debug = False
        if rank==0:
            self.algorithm_initialize()
            logger.info('Beginning main algorithm loop')
        #Begin main algorithm loop
        self.convergence = False
        convergence=False
        while not convergence:
            if rank==0:
                pop = self.population
                offspring = self.generation_set(pop)
                # Identify the individuals with an invalid fitness
                invalid_ind = [ind for ind in offspring if ind.fitness==0]
                #Evaluate the individuals with invalid fitness
                self.output.write('\n--Evaluate Structures--\n')
                ntimes=int(math.ceil(float(len(invalid_ind))/float(comm.Get_size())))
                nadd=int(ntimes*comm.Get_size()-len(invalid_ind))
                maplist=[[] for n in range(ntimes)]
                strt=0
                for i in range(len(maplist)):
                    maplist[i]=[[self,indi] for indi in invalid_ind[strt:comm.Get_size()+strt]]
                    strt+=comm.Get_size()
                for i in range(nadd):
                    maplist[len(maplist)-1].append([None,None])
            else:
                ntimes=None
            ntimes = comm.bcast(ntimes,root=0)
            outs=[]
            for i in range(ntimes):
                if rank==0:
                    one=maplist[i]
                else:
                    one=None
                ind =comm.scatter(one,root=0)
                out = switches.fitness_switch(ind)
                outt = comm.gather(out,root=0)
                if rank==0:
                    outs.extend(outt)
            if rank==0:
                for i in range(len(invalid_ind)):
                    invalid_ind[i] = outs[i][0]
                for i in range(len(outs)):
                    self.output.write(outs[i][1])
                pop.extend(invalid_ind)
                pop = self.generation_eval(pop)
                self.write()
            convergence =comm.bcast(self.convergence, root=0)
        
        if rank==0:
            logger.info('Run algorithm stats')
            end_signal = self.algorithm_stats(self.population)
        else:
            end_signal = None
        end_signal = comm.bcast(end_signal, root=0)
        return end_signal    
    
    def algorithm_par_mp1(self):
        """Subprogram for running parallel version of GA
        Requires MPI4PY"""
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        if rank==0:
            if 'MA' in self.debug:
                debug = True
            else:
                debug = False
            self.algorithm_initialize()
        self.convergence = False
        convergence = False
        while not convergence:
            if rank==0:
                self.calc = tools.setup_calculator(self) #Set up calculator for atomic structures
                #Set up calculator for fixed region calculations
                if self.fixed_region:
                    self.static_calc = self.calc #May need to copy this
                    self.calc = tools.setup_fixed_region_calculator(self)
                pop = self.population
                offspring = self.generation_set(self,pop)
                # Identify the individuals with an invalid fitness
                invalid_ind = [ind for ind in offspring if ind.energy==0]
                #Evaluate the individuals with invalid fitness
                self.output.write('\n--Evaluate Structures--\n')
                proc_dist = int(comm.Get_size()/self.n_proc_eval)
                ntimes=int(math.ceil(float(len(invalid_ind))/float(proc_dist)))
                nadd=int(ntimes*proc_dist-len(invalid_ind))
                maplist=[[] for n in range(ntimes)]
                strt=0
                for i in range(len(maplist)):
                    maplist[i]=[[self,indi] for indi in invalid_ind[strt:proc_dist+strt]]
                    strt+=proc_dist
                for i in range(nadd):
                    maplist[len(maplist)-1].append([None,None])
                masterlist = [i*self.n_proc_eval for i in range(proc_dist)]
            else:
                ntimes=None
                masterlist = None
            ntimes = comm.bcast(ntimes,root=0)
            outs=[]
            for i in range(ntimes):
                if rank==0:
                    one=maplist[i]
                    for j in range(len(one)):
                        comm.send(ind, dest=1, tag=11)
                elif rank == 1:
                    ind = comm.recv(source=0, tag=11)
                else:
                    one=None
                ind =comm.scatter(one,root=0)
                out = switches.fitness_switch(ind)
            else:
                invalid_ind=[]
            poorlist = []
            for i in range(len(invalid_ind)):
                if self.fitness_scheme=='STEM_Cost':
                    if self.stem_coeff==None:
                        ind = invalid_ind.pop()
                        from MAST.structopt.tools.StemCalc import find_stem_coeff
                        outs = find_stem_coeff(self,ind)
                        ind = outs[1]
                        self.stem_coeff = outs[0]
                        self.output.write('stem_coeff Calculated to be: '+repr(self.stem_coeff)+'\n')
                        pop.append(ind)
                ind=invalid_ind[i]
                if 'MA' in self.debug: write_xyz(self.debugfile,ind[0],'Individual to fitness_switch')
                outs = switches.fitness_switch([self,ind])
                self.output.write(outs[1])
                invalid_ind[i]=outs[0]
                if invalid_ind[i].energy == float('inf'):
                    poorlist.append(i)
                    self.output.write('Removing infinite energy individual '+repr(ind.history_index)+'\n')
                elif invalid_ind[i].energy == float('nan'):
                    poorlist.append(i)
                    self.output.write('Removing nan energy individual '+repr(ind.history_index)+'\n')
                self.output.flush()
            if len(poorlist) != 0:
                poorlist.sort(reverse=True)
                for one in poorlist:
                    del invalid_ind[one]
            if rank==0:
                pop.extend(invalid_ind)
                pop = self.generation_eval(pop)
            convergence = comm.bcast(self.convergence, root=0)
        
        if rank==0:
            end_signal = self.algorithm_stats(self.population)
        else:
            end_signal = None
        end_signal = comm.bcast(end_signal, root=0)
        return end_signal    
    
    def algorithm_par_mp(self):
        """Subprogram for running parallel version of GA
        Requires MPI4PY"""
        global logger
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        if rank==0:
            if 'MA' in self.debug:
                debug = True
            else:
                debug = False
            self.algorithm_initialize()
        self.convergence = False
        convergence = False
        while not convergence:
            if rank==0:
                pop = self.population
                offspring = self.generation_set(self,pop)
                # Identify the individuals with an invalid fitness
                invalid_ind = [ind for ind in offspring if ind.energy==0]
                #Evaluate the individuals with invalid fitness
                self.output.write('\n--Evaluate Structures--\n')
            else:
                invalid_ind=[]
            for i in range(len(invalid_ind)):
                if self.fitness_scheme=='STEM_Cost':
                    if self.stem_coeff==None:
                        ind = invalid_ind.pop()
                        from MAST.structopt.tools.StemCalc import find_stem_coeff
                        outs = find_stem_coeff(self,ind)
                        ind = outs[1]
                        self.stem_coeff = outs[0]
                        self.output.write('stem_coeff Calculated to be: '+repr(self.stem_coeff)+'\n')
                        pop.append(ind)
                ind=invalid_ind[i]
                if 'MA' in self.debug: write_xyz(self.debugfile,ind[0],'Individual to fitness_switch')
                outs = switches.fitness_switch([self,ind])
                self.output.write(outs[1])
                invalid_ind[i]=outs[0]
                self.output.flush()
            if rank==0:
                pop.extend(invalid_ind)
                pop = self.generation_eval(pop)
                self.write()
            convergence = comm.bcast(self.convergence, root=0)
        
        if rank==0:
            end_signal = self.algorithm_stats(self.population)
        else:
            end_signal = None
        end_signal = comm.bcast(end_signal, root=0)
        return end_signal    
    
    def algorithm_serial(self):
        """Subprogram to run the optimizer in serial"""
        global logger
        self.algorithm_initialize()
        logger.info('Beginning main algorithm loop')
        #Begin main algorithm loop
        while not self.convergence:
            logger.info('Setup calculator for generation {0}'.format(self.generation))
            pop = self.population
            logger.info('Setup population')
            offspring = self.generation_set(pop)
            # Identify the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if ind.energy==0]
            #DEBUG: Write first invalid ind
            if 'MA' in self.debug:
                logger.info('Identified {0} structures with energy=0'.format(len(invalid_ind)))
                inp_out.write_xyz(self.debugfile,invalid_ind[0][0],\
                'First from Invalid_ind list '+repr(invalid_ind[0].energy))
                #DEBUG: Write first invalid ind in solid
                if self.structure=='Defect' or self.structure=='Surface':
                    sols = invalid_ind[0][0].copy()
                    sols.extend(invalid_ind[0].bulki)
                    inp_out.write_xyz(self.debugfile,sols,'First from Invalid-ind + Bulki')
                    sols = invalid_ind[0][0].copy()
                    sols.extend(invalid_ind[0].bulko)
                    inp_out.write_xyz(self.debugfile,sols,'First from Invalid-ind + Bulko')
            self.output.write('\n--Evaluate Structures--\n')
            logger.info('Evaluating fitness of structures')
            for i in range(len(invalid_ind)):
                ind = invalid_ind[i]
                if 'MA' in self.debug:
                    write_xyz(self.debugfile,ind[0],'Individual to fitness_switch')
                outs = switches.fitness_switch([self,ind])
                self.output.write(outs[1])
                invalid_ind[i] = outs[0]
            pop.extend(invalid_ind)
            logger.info('Eval Population')
            pop = self.generation_eval(pop)
            self.write()
            
        logger.info('Run algorithm stats')
        end_signal = self.algorithm_stats(self.population)
        return end_signal    
    
    def algorithm_stats(self,pop):
        self.output.write('\n----- Algorithm Stats -----\n')
        cxattempts = 0
        cxsuccess = 0
        for ats,sus in self.CXs:
            cxattempts+=ats
            cxsuccess+=sus
        mutslist = [[0,0] for one in self.mutation_options]
        for one in self.Muts:
            for i in range(len(mutslist)):
                mutslist[i][0]+=one[i][0]
                mutslist[i][1]+=one[i][1]
        self.output.write('Total Number of Evaluations : '+repr(sum(self.Evaluations))+'\n')
        self.output.write('Average Number of Evaluations per Generation : '+repr(
            float(sum(self.Evaluations))/float(len(self.Evaluations)))+'\n')
        ttseconds = self.Runtimes[-1]-self.Runtimes[0]
        deltats = []
        for i in range(1,len(self.Runtimes)):
            deltats.append(self.Runtimes[i]-self.Runtimes[i-1])
        maxs = [(ind,value) for ind,value in enumerate(deltats) if value==max(deltats)]
        maxtgen=maxs[0][0]-1
        maxts=maxs[0][1]
        mins = [(ind,value) for ind,value in enumerate(deltats) if value==min(deltats)]
        mintgen=mins[0][0]-1
        mints=mins[0][1]
        avgts = sum(deltats)/len(deltats)
        self.output.write('Total Length of GA run : '+tools.convert_time(ttseconds)+'\n')
        self.output.write('Average time per generation : '+tools.convert_time(avgts)+'\n')
        self.output.write('Maximum time for generation '+repr(maxtgen)+' : '+tools.convert_time(maxts)+'\n')
        self.output.write('Minimum time for generation '+repr(mintgen)+' : '+tools.convert_time(mints)+'\n')
        self.output.write('Attempted Crossovers: ' + repr(cxattempts)+'\n')
        self.output.write('Successful Crossovers: ' + repr(cxsuccess)+'\n')
        self.output.write('Mutations:\n')
        i=0
        for opt in self.mutation_options:
            self.output.write('    Attempted ' + opt + ' : ' + repr(mutslist[i][0]) + '\n')
            self.output.write('    Successful ' + opt + ' : ' + repr(mutslist[i][1]) + '\n')
            i+=1
        if self.best_inds_list:
            BESTS=tools.BestInds(pop,self.BESTS,self,writefile=True)
        self.close_output()
        end_signal='Genetic Algorithm Finished'
        return end_signal
    
    def check_pop(self, pop):
        # Gather all the energies/fitnesses
        if self.output_format=='totalenergy':
            complist = [ind.energy for ind in pop]
        elif self.output_format=='formation_energy':
            complist = []
            for ind in pop:
                solid=Atoms()
                solid.extend(ind[0])
                solid.extend(ind.bulki)
                energy=ind.energy
                for sym,c,m,u in self.atomlist:
                    nc=len([atm for atm in solid if atm.symbol==sym])
                    energy-= float(nc)*float(u)
                complist.append(energy)
        elif self.output_format=='formation_energy_per_int':
            complist = []
            for ind in pop:
                solid=Atoms()
                solid.extend(ind[0])
                solid.extend(ind.bulki)
                energy=ind.energy
                for sym,c,m,u in self.atomlist:
                    nc=len([atm for atm in solid if atm.symbol==sym])
                    energy-= float(nc)*float(u)
                energy=energy/self.natoms
                complist.append(energy)
        elif self.output_format=='formation_energy2':
            complist = [(ind.energy - ind.purebulkenpa*(ind[0].get_number_of_atoms()+ind.bulki.get_number_of_atoms()))/(ind[0].get_number_of_atoms()+ind.bulki.get_number_of_atoms()-ind.natomsbulk) for ind in pop]
        elif self.output_format=='energy_per_atom':
            complist = [ind.energy/(ind[0].get_number_of_atoms()+ind.bulki.get_number_of_atoms()) for ind in pop]
        else:
            complist = [ind.fitness for ind in pop]
    
        # Calcluate and print the Stats
        length = len(pop)
        mean = sum(complist) / length
        complist.sort()
        medium = complist[length/2]
        sum2 = sum(x*x for x in complist)
        std = abs(sum2 / length - mean**2)**0.5
        mine = min(complist)
        maxe = max(complist)
    
        self.output.write('\n----Stats----\n')
        self.output.write('  Min '+repr(mine)+'\n')
        self.output.write('  Max '+repr(maxe)+'\n')
        self.output.write('  Avg '+repr(mean)+'\n')
        self.output.write('  Medium '+repr(medium)+'\n')
        self.output.write('  Std '+repr(std)+'\n')
        self.output.write('  Genrep '+repr(self.genrep)+'\n')
        self.summary.write('{Gen} {Fitmin} {Fitavg} {Fitmed} {Fitmax} {Std} {time}\n'.format(
            Gen=str(self.generation),Fitmin=repr(mine),Fitavg=repr(mean),Fitmed=repr(medium),Fitmax=repr(maxe),
            Std=repr(std),time=repr(time.asctime( time.localtime(time.time()) ))))
    
        # Set new index values and write population
        index1 = 0
        for ind in pop:
            ind.index=index1
            index1+=1
        inp_out.write_pop(self,pop)
    
        if self.allenergyfile:
            for ind in pop:
                self.tenergyfile.write(repr(ind.energy)+' ')
            self.tenergyfile.write('\n')
    
        #Check Convergence of population based on fitness
        fitnesses = [ind.fitness for ind in pop]
        popmin = min(fitnesses)
        popmean = sum(fitnesses) / len(fitnesses)
        sum2 = sum(x*x for x in fitnesses)
        popstd = abs(sum2 / len(fitnesses) - popmean**2)**0.5
        convergence = False
        if self.convergence_scheme=='gen_rep_min':
            if self.generation < self.maxgen:
                if self.generation == 0:
                    self.minfit = popmin
                if abs(self.minfit - popmin) < self.tolerance:
                    self.genrep += 1
                else:
                    self.genrep = 0
                    self.minfit = popmin
                if self.genrep > self.reqrep:
                    self.convergence = True
            else:
                self.convergence = True
        elif self.convergence_scheme=='gen_rep_avg':
            if self.generation < self.maxgen:
                if self.generation == 0:
                    self.minfit = popmean
                if abs(self.minfit - popmean) < self.tolerance:
                    self.genrep += 1
                else:
                    self.genrep = 0
                    self.minfit = popmean
                if self.genrep > self.reqrep:
                    self.convergence = True
            else:
                self.convergence = True
        elif self.convergence_scheme=='std':
            if popstd < self.tolerance:
                self.convergence = True
        else:#self.convergence_scheme=='max_gen':
            if self.generation > self.maxgen:
                self.convergence = True
    
        #Flush output to files
        self.output.flush()
        self.summary.flush()
        for ind in self.files:
            ind.flush()
        if self.indiv_defect_write:
            for ind in self.ifiles:
                ind.flush()
        if self.genealogy: self.Genealogyfile.flush()
        if self.allenergyfile: self.tenergyfile.flush()
        if 'MA' in self.debug: self.debugfile.flush()
        if self.fingerprinting: 
            self.fpfile.flush()
            self.fpminfile.flush()
        return self

    def close_output(self):
        localtime = time.asctime( time.localtime(time.time()) )
        for ind in self.files:
            ind.close()
        self.output.write('Local time : '+repr(localtime)+'\n')
        self.output.write('End of Execution\n')
        self.summary.close()
        self.output.close()
        if self.genealogy: self.Genealogyfile.close()
        if self.allenergyfile: self.tenergyfile.close()
        if 'MA' in self.debug: self.debugfile.close()
        if self.fingerprinting: 
            self.fpfile.close()
            self.fpminfile.close()
    
    def generation_eval(self, pop):
        global logger
        emx = max(ind.energy for ind in pop)
        emn = min(ind.energy for ind in pop)
        for ind in pop:
            ind.tenergymx = emx
            ind.tenergymin = emn
        #DEBUG: Write relaxed individual
        if 'MA' in self.debug:
            if self.generation > 0: 
                inp_out.write_xyz(self.debugfile,pop[self.nindiv][0],\
                'First Relaxed Offspring '+repr(pop[self.nindiv-1].energy))    
                #DEBUG: Write relaxed ind in solid
                if self.structure=='Defect' or self.structure=='Surface':
                    inp_out.write_xyz(self.debugfile,pop[self.nindiv].bulki,\
                    'First Relaxed bulki '+repr(pop[self.nindiv-1].energy))    
                    sols = pop[self.nindiv][0].copy()
                    sols.extend(pop[self.nindiv].bulki)
                    inp_out.write_xyz(self.debugfile,sols,'First from Invalid-ind + Bulki '+\
                    repr(pop[self.nindiv].energy))
                    sols = pop[self.nindiv][0].copy()
                    sols.extend(pop[self.nindiv].bulko)
                    inp_out.write_xyz(self.debugfile,sols,\
                    'First from Invalid-ind + Bulko '+repr(pop[self.nindiv].energy))
        if self.generation==0:
            logger.info('Initializing Bests list')
            self.BESTS = list()
        if self.best_inds_list:
            self.BESTS = tools.BestInds(pop,self.BESTS,self,writefile=True)
        # Determine survival based on fitness predator
        if 'lambda,mu' not in self.algorithm_type:
            pop = tools.get_best(pop, len(pop))
        if self.fingerprinting:
            logger.info('Writing fingerprint files')
            for one in pop:
                self.fpfile.write(repr(fingerprinting.fingerprint_dist(
                    pop[0].fingerprint,one.fingerprint))+' '+repr(one.energy)+' ')
            self.fpfile.write('\n')
            self.fpminfile.write(repr(pop[0].fingerprint)+'\n')
            self.fpminfile.write(repr(pop[0].energy)+'\n')
        nevals = len(pop)/2
        if self.generation !=0:
            logger.info('Applying predator')
            pop = predator_switch(pop,self)
        else:
            self.genrep = 0
            self.minfit = 0
        # Evaluate population
        logger.info('Checking population for convergence')
        self.check_pop(pop)
        #Update general output tracking
        if self.generation !=0:
            histlist = []
            for ind in pop:
                histlist.append(ind.history_index)
            self.Runtimes.append(time.time())
            self.Evaluations.append(nevals)
            cxsuccess = 0
            mutsuccess = []
            for one in histlist:
                if '+' in one:
                    cxsuccess +=1
                if 'm' in one:
                    mutsuccess.append(one)
            self.CXs.append((self.cxattempts,cxsuccess))
            mutslist = [[0,0] for one in self.mutation_options]
            for one in mutsuccess:
                for two, opt in self.mutattempts:
                    if one==two:
                        index = [ind for ind,value in enumerate(
                            self.mutation_options) if value==opt][0]
                        mutslist[index][1]+=1
            for one,opt in self.mutattempts:
                index = [ind for ind,value in enumerate(
                self.mutation_options) if value==opt][0]
                mutslist[index][0]+=1
            self.Muts.append(mutslist)
            self.output.write('\n----- Generation Stats -----\n')
            self.output.write('Attempted Crossovers: ' + repr(self.cxattempts)+'\n')
            self.output.write('Successful Crossovers: ' + repr(cxsuccess)+'\n')
            self.output.write('Mutations:\n')
            i=0
            for opt in self.mutation_options:
                self.output.write('    Attempted ' + opt + ' : ' + repr(mutslist[i][0]) + '\n')
                self.output.write('    Successful ' + opt + ' : ' + repr(mutslist[i][1]) + '\n')
                i+=1
        self.generation += 1
        # Set new index values
        index1 = 0
        for ind in pop:
            ind.index = index1
            index1+=1
        index1 = 0
        if not self.convergence:
            try:
                self.calc.clean()
                if self.fixed_region:
                    self.static_calc.clean()
            except:
                pass
            if self.lammps_keep_files:
                try:
                    if self.parallel and ('Island_Method' not in self.algorithm_type):
                        nproc = MPI.COMM_WORLD.Get_size()
                        lmpfilepath = os.path.dirname(self.calc.tmp_dir)
                        for proc in range(nproc):
                            pth = os.path.join(lmpfilepath,'rank-{0}'.format(proc))
                            fls = [fl for fl in os.listdir(pth) if '.' not in fl]
                            for one in fls:
                                os.remove(pth+'/'+one)
                    else:
                        fls = [fl for fl in os.listdir(self.calc.tmp_dir) if '.' not in fl]
                        for one in fls:
                            os.remove(self.calc.tmp_dir+'/'+one)
                except Exception, e:
                    logger.error('Print issue in removing files {0}.'.format(e),exc_info=True)
                    print str(e)
                    pass
        self.population = pop
        return pop
    
    def generation_set(self,pop):
        global logger
        self.calc = tools.setup_calculator(self) #Set up calculator for atomic structures
        #Set up calculator for fixed region calculations
        if self.fixed_region:
            self.static_calc = self.calc #May need to copy this
            self.calc = tools.setup_fixed_region_calculator(self)
        self.output.write('\n-------- Generation '+repr(self.generation)+' --------\n')
        self.files[self.nindiv].write('Generation '+str(self.generation)+'\n')
        if len(pop) == 0:
            logger.info('Initializing structures')
            offspring = self.initialize_structures()
            self.population = offspring
        else:
            for i in range(len(pop)):
                # Reset History index
                pop[i].history_index=repr(pop[i].index)
            # Select the next generation individuals
            offspring = switches.selection_switch(pop, self.nindiv,
                        self.selection_scheme, self)
            # Clone the selected individuals
            offspring=[off1.duplicate() for off1 in offspring]
            # Apply crossover to the offspring
            self.output.write('\n--Applying Crossover--\n')
            cxattempts = 0
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cxpb:
                    child1,child2 = switches.crossover_switch(child1, child2, self)
                    cxattempts+=2
            self.cxattempts=cxattempts
            #DEBUG: Write first child
            if 'MA' in self.debug: 
                inp_out.write_xyz(self.debugfile,offspring[0][0],'First Child '+
                    repr(offspring[0].history_index))
            # Apply mutation to the offspring
            self.output.write('\n--Applying Mutation--\n')
            mutattempts = []
            muts = []
            for mutant in offspring:
                if random.random() < self.mutpb:
                    if self.mutant_add:
                        mutant = mutant.duplicate()
                    mutant, optsel = switches.moves_switch(mutant,self)
                    mutattempts.append([mutant.history_index,optsel])
                    if self.mutant_add:
                        muts.append(mutant)
            if self.mutant_add:
                offspring.extend(muts)
            self.mutattempts=mutattempts
            #DEBUG: Write first offspring
            if 'MA' in self.debug: 
                inp_out.write_xyz(self.debugfile,muts[0][0],'First Mutant '+\
                repr(muts[0].history_index))
        if 'stem' in self.fitness_scheme:
            if self.stem_coeff==None:
                logger.info('Setting STEM coeff (alpha)')
                ind = offspring.pop()
                from MAST.structopt.tools.StemCalc import find_stem_coeff
                outs = find_stem_coeff(self,ind)
                ind = outs[1]
                ind.fitness = 0
                ind.energy = 0
                self.stem_coeff = outs[0]
                logger.info('STEM Coeff = {0}'.format(self.stem_coeff))
                self.output.write('stem_coeff Calculated to be: '+repr(self.stem_coeff)+'\n')
                offspring.append(ind)
                  
        return offspring
    
    def initialize_structures(self):
        global logger
        self.output.write('\n----Initialize Structures----\n')
        #self.logger.info('Initializing Structures')
        # Initialize Population - Generate a list of ncluster individuals
        # Set bulk and index atributes
        if self.restart:
            logger.info('Loading previous population')
            pop = generate.get_restart_population(self)
        else:
            logger.info('Generating new population')
            pop = generate.get_population(self)
        if 'MA' in self.debug: 
            inp_out.write_xyz(self.debugfile,pop[0][0],'First Generated Individual')
        #Use if concentration of interstitials is unknown
        if self.swaplist:
            mutlist=self.mutation_options
            self.mutation_options=['IntSwapLocal']
            for i in range(len(pop)):
                one = pop[i]
                one.swaplist=copy.deepcopy(swaplist)
                if random.random() < 0.25:
                    pop[i], opt = switches.moves_switch(one,self)
            self.mutation_options=mutlist
        # Write Initial structures to files
        self.output.write('\n---Starting Structures---\n')
        inp_out.write_pop(self, pop)
        #Print number of atoms to Summary file
        if self.structure=='Defect' or self.structure=='Surface':
            natstart=len(pop[0][0])+len(pop[0].bulki)
        else:
            natstart=len(pop[0][0])
        if not self.restart:
            if self.structure=='Defect':
                self.summary.write('Defect Run Pure Bulk Energy per Atom : '+repr(self.purebulkenpa)+'\n')
            else:
                self.summary.write(self.structure+' Run : '+repr(0)+'\n')
            self.summary.write('Natoms '+repr(natstart)+ '\n')
            #Print data headers to summary file
            self.summary.write('Generation Fitmin Fitavg Fitmedium Fitmax std time \n')
        offspring=pop
        #pop=[]
        #BESTS=[]
        return offspring
    
    def run(self):
        global logger
        cwd=os.getcwd()
        try:
            if self.parallel:
                MPI.COMM_WORLD.Set_errhandler(MPI.ERRORS_ARE_FATAL) 
                if 'Island_Method' in self.algorithm_type:
                    logger.info('Running Island Method algorithm')
                    self.algorithm_island()
                    if self.postprocessing:
                        logger.info('Running Post-processing')
                        rank = MPI.COMM_WORLD.Get_rank()
                        path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                        os.chdir(path)
                        if self.genealogytree:
                            pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms,loggername=self.loggername)
                        else:
                            pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms, loggername=self.loggername)
                        os.chdir(cwd)
                    if self.lattice_concentration:
                        if self.structure=='Defect':
                            logger.info('Running lattice concentration check')
                            rank = MPI.COMM_WORLD.Get_rank()
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                            os.chdir(path)
                            if self.best_inds_list:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+self.filename+'.xyz'))
                            else:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                            os.chdir(cwd)
                elif 'MP_Eval' in self.algorithm_type:
                    logger.info('Running multiprocessor parallel algorithm')
                    self.algorithm_par_mp()
                    rank = MPI.COMM_WORLD.Get_rank()
                    if rank==0:
                        if self.postprocessing:
                            logger.info('Running Post-processing')
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                            os.chdir(path)
                            if self.genealogytree:
                                pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms,loggername=self.loggername)
                            else:
                                pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms, loggername=self.loggername)
                            os.chdir(cwd)
                    if self.lattice_concentration:
                        if self.structure=='Defect':
                            logger.info('Running lattice concentration check')
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                            os.chdir(path)
                            if self.best_inds_list:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+self.filename+'.xyz'))
                            else:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                            os.chdir(cwd)
                else:
                    logger.info('Running parallel algorithm')
                    self.algorithm_parallel()
                    rank = MPI.COMM_WORLD.Get_rank()
                    if rank==0:
                        if self.postprocessing:
                            logger.info('Running Post-processing')
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                            os.chdir(path)
                            if self.genealogytree:
                                pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms)
                            else:
                                pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms)
                            os.chdir(cwd)
                        if self.lattice_concentration:
                            if self.structure=='Defect':
                                logger.info('Running lattice concentration check')
                                path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                                os.chdir(path)
                                if self.best_inds_list:
                                    pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+self.filename+'.xyz'))
                                else:
                                    pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                                os.chdir(cwd)
            else:
                logger.info('Running serial algorithm')
                self.algorithm_serial()
                rank = 0
                if self.postprocessing:
                    logger.info('Running Post-processing')
                    path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                    os.chdir(path)
                    if self.genealogytree:
                        pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms)
                    else:
                        pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms)
                    os.chdir(cwd)
                if self.lattice_concentration:
                    if self.structure=='Defect':
                        logger.info('Running lattice concentration check')
                        path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(self.filename,rank))
                        os.chdir(path)
                        if self.best_inds_list:
                            pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+self.filename+'.xyz'))
                        else:
                            pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                        os.chdir(cwd)
        except Exception, e:
            logger.error('Error in execution: {0}'.format(e),exc_info=True)
            print '********ERROR IN EXECUTION********'
            print str(e)
            print 'CLEANING UP FILES'
            try:
                self.close_output()
            except:
                pass
            print 'EXITING PROGRAM'
    
    def read(self,optfile):
        parameters = inp_out.read_parameter_input(optfile,True)
        self.__dict__.update(parameters)
        outdict = inp_out.restart_output(self)
        self.__dict__.update(outdict)
        poplist = []
        for indfile in self.population:
            ind = inp_out.read_individual(indfile)
            poplist.append(ind)
        self.population = poplist
        bestlist = []
        for bestfile in self.BESTS:
            ind = inp_out.read_individual(bestfile)
            bestlist.append(ind)
        self.BESTS = bestlist
        self.restart = True
        return self
        
    def write(self,filename=None, restart=True):
        if filename:
            inp_out.write_optimizer(self, filename, restart)
        else:
            inp_out.write_optimizer(self, self.optimizerfile, restart)
        return
    
if __name__ == "__main__":
    import sys
    input = sys.argv[1]
    A=Optimizer(input)
    A.run()
