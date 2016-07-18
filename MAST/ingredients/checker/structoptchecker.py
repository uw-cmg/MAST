import os
import shutil
import sys
import time
import logging
import copy
import subprocess
from MAST.utility import loggerutils
from MAST.utility import MASTFile
from MAST.utility import dirutil
from MAST.utility import MASTError
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import LammpsChecker
from MAST.submit import queue_commands
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from MAST.structopt.Optimizer import Optimizer
from MAST.structopt import post_processing as pp
from MAST.structopt import inp_out
from MAST.structopt.switches import fitness_switch

class StructoptChecker(BaseChecker):
    
    def __init__(self, **kwargs):
        """Please not modify this init method."""
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)
        #print program_keys
        self.structopt_parameters = self._structopt_get_non_mast_keywords()
        self.potential_file = None
        if 'potential_file' in self.structopt_parameters:
            if self.structopt_parameters['potential_file'] != 'None':
                if os.path.isfile(self.structopt_parameters['potential_file']):
                    self.potential_file = self.structopt_parameters['potential_file']
                else:
                    raise MASTError(self.__class__.__name__,
                        "Cannot find specified potential file for LAMMPS: %s" % self.structopt_parameters['potential_file'])
        if "MAST" in self.structopt_parameters['calc_method']:
            if "exec_mast" not in self.structopt_parameters:
                raise MASTError(self.__class__.__name__,
                    "Must specify command for Mast execution in structopt. Use exec_mast command")

    def _structopt_get_non_mast_keywords(self):
        """Get the StructOpt keywords and make a dictionary."""
        input_dict=dict()
        allowedpath = os.path.join(dirutil.get_mast_install_path(),
                        'ingredients', 
                        'programkeys','structopt_allowed_keywords.py')
        allowed_list = self._structopt_get_allowed_keywords(allowedpath)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                keytry = key
                if not (keytry in allowed_list):
                    self.logger.warning("Ignoring program key %s for INPUT. To allow this keyword, add it to %s" % (keytry, allowedpath))
                else:
                    try:
                        #Convert numbers to floats
                        value = float(value)
                    except:
                        try:
                            #Convert list, tuple, and booleans
                            if value != 'calc_method':
                                value = eval(value)
                            else:
                                value = value.strip()
                        except:
                            try:
                                #Leave remaining as strings
                                value = value.strip()
                            except:
                                print 'Trouble with input line: ', value
                    input_dict[keytry]=value
        if 'structopt_input_file' in input_dict:
            self.logger.info('Requested read from StructOpt input file at : {0}'.format(input_dict['structopt_input_file']))
            new_dict = inp_out.read_parameter_input(input_dict['structopt_input_file'],False)
            #Overwrite parameters from input file with parameters in mast input
            for key, value in input_dict.iteritems():
                new_dict[key] = value
            input_dict = new_dict
            self.logger.info('New input_dict as read from StructOpt input file : {0}'.format(input_dict))
        # Check for potential file
        if 'potential_file' in input_dict:
            if 'pot_file' not in input_dict:
                input_dict['pot_file']
        return input_dict

    def _structopt_get_allowed_keywords(self, allowedpath):
        """Get allowed vasp keywords.
            Args:
                allowedpath <str>: file path for allowed lammps keywords
        """
        allowed = MASTFile(allowedpath)
        allowed_list=list()
        for entry in allowed.data:
            allowed_list.append(entry.strip())
        return allowed_list
    
    def transfer_restart_files(self, input_dict):
        if 'filename' in input_dict:
            fname = '{0}-rank0'.format(input_dict['filename'])
        else:
            fname = 'Output-rank0'
        if ('restart' in input_dict) and (input_dict['restart']):
            try:
                restartdir = os.path.join(self.keywords['name'],fname)
                files = os.listdir(restartdir)
            except:
                raise MASTError(self.__class__.__name__,
                    "Cannot find directory for restart: {0}".format(restartdir))
            return input_dict
        if input_dict['restart_optimizer']:
            if 'ARCHIVE' in input_dict['structopt_input_file']:
                replaceflag = True
            else:
                replaceflag = False
            path = os.path.join(self.keywords['name'],fname)
            if not os.path.exists(path):
                os.mkdir(path)
            outputoptions = ['optimizerfile', 'tenergyfile', 'fpfile',
                            'fpminfile', 'debugfile', 'Genealogyfile', 'summary']
            for one in outputoptions:
                if one in input_dict:
                    if (input_dict[one] != None):
                        bname = os.path.basename(input_dict[one])
                        npath = os.path.join(path,bname)
                        if replaceflag:
                            input_dict[one] = input_dict[one].replace('SCRATCH','ARCHIVE')
                        shutil.copyfile(input_dict[one], npath)
                        input_dict[one] = npath
            flist = list()
            for one in input_dict['files']:
                bname = os.path.basename(one)
                npath = os.path.join(path,bname)
                if replaceflag:
                    one = one.replace('SCRATCH','ARCHIVE')
                shutil.copyfile(one, npath)
                flist.append(npath)
            input_dict['files'] = flist
            iflist = list()
            if input_dict['ifiles']:
                for one in input_dict['ifiles']:
                    bname = os.path.basename(one)
                    npath = os.path.join(path,bname)
                    if replaceflag:
                        one = one.replace('SCRATCH','ARCHIVE')
                    shutil.copyfile(one, npath)
                    iflist.append(npath)
            input_dict['ifiles'] = iflist
            bname = os.path.basename(input_dict['output'])
            npath = os.path.join(self.keywords['name'],bname)
            if replaceflag:
                input_dict['output'] = input_dict['output'].replace('SCRATCH','ARCHIVE')
            fout = open(input_dict['output'],'r')
            lines = fout.readlines()
            fout.close()
            nlines = []
            for line in lines:
                if "End of Execution" not in line:
                    nlines.append(line)
            foutn = open(npath,'a')
            for line in nlines:
                foutn.write(line)
            foutn.close()
            #shutil.copyfile(input_dict['output'], npath)
            input_dict['output'] = npath
            path = os.path.join(path, 'Restart-files')
            if not os.path.exists(path):
                os.mkdir(path)
            popfiles = list()
            for one in input_dict['population']:
                bname = os.path.basename(one)
                npath = os.path.join(path,bname)
                if replaceflag:
                    one = one.replace('SCRATCH','ARCHIVE')
                shutil.copy(one, npath)
                popfiles.append(npath)
            input_dict['population'] = popfiles
            bestfiles = list()
            for one in input_dict['BESTS']:
                bname = os.path.basename(one)
                npath = os.path.join(path,bname)
                if replaceflag:
                    one = one.replace('SCRATCH','ARCHIVE')
                shutil.copy(one, npath)
                bestfiles.append(npath)
            input_dict['BESTS'] = bestfiles
            return input_dict
    
    def is_frozen(self):
        """Check if StructOpt calculation is frozen.
        """
        if 'filename' in self.structopt_parameters:
            outfile = self.structopt_parameters['filename']
        else:
            outfile = "Output.txt"
        if not os.path.isfile(os.path.join(self.keywords['name'], outfile)):
            return False
        return BaseChecker.is_frozen(self,outfile)
    
    def is_complete(self):
        if "MAST" in self.structopt_parameters['calc_method']:
            if "LOOPED" in self.structopt_parameters['calc_method']:
                return self.is_complete_looped()
            else:
                return self.is_complete_one_gen()
        else:
            return self.is_complete_whole()
        
    def is_complete_whole(self):
        """Check if StructOpt calculation is complete.
        """
        if 'filename' in self.structopt_parameters:
            fname = '{0}-rank0.txt'.format(self.structopt_parameters['filename'])
            opath = os.path.join(self.keywords['name'],fname)
        else:
            opath = os.path.join(self.keywords['name'],'Output-rank0.txt')
        if not os.path.isfile(opath):
            self.logger.info("No StructOpt output at %s; not complete." % opath)
            return False
        reachgrep=subprocess.Popen('grep "End of Execution" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            self.logger.info("StructOpt output file at %s shows calculation is incomplete." % opath)
            return False
        else:
            self.logger.info("StructOpt output file at %s shows calculation has finished." % opath)
            return True
    
    def is_complete_one_gen(self):
        """Check if StructOpt single generation is complete.
        """
        if 'filename' in self.structopt_parameters:
            opath = os.path.join(self.keywords['name'],self.structopt_parameters['filename'])
        else:
            opath = os.path.join(self.keywords['name'],'Output-rank0.txt')
#         if not os.path.isfile(opath):
#             self.logger.info("No StructOpt output at %s; not complete." % opath)
#             return False
        reachgrep=subprocess.Popen('grep "End of Execution" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            reachedend=False
        else:
            self.logger.info("StructOpt output file at %s shows calculation has finished." % opath)
            return True
        if 'filename' in self.structopt_parameters:
            optifile = os.path.join(os.path.join(self.keywords['name'], '{0}-rank0'.format(self.structopt_parameters['filename']),'Optimizer-restart-file.txt'))
        else:
            optifile = os.path.join(os.path.join(self.keywords['name'], 'Output-rank0','Optimizer-restart-file.txt'))
        ingpath = self.keywords['name']
        if not os.path.isfile(optifile):
            os.chdir(ingpath)
            inputfile = os.path.join(ingpath,"structoptinput.txt")
            Opti = Optimizer(inputfile)
            Opti.algorithm_initialize()
            self.logger.info("Opti initialized")
            self.logger.info('Creating files for initial population')
        else:
            if not self.subfolders_complete():
                self.logger.info("Subfolders not complete. No further checks.")
                self.flip_status_back()
                return False
            else: 
                self.logger.info("Evaluating completed structures.")
                Opti = Optimizer(input=optifile)
                self.logger.info('Optimizer read optimizer restart file')
                self.logger.info("Generation %i" % Opti.generation)
                Opti = self.evaluate_structures_for_fitness(Opti)
                if Opti.convergence:
                    os.chdir(ingpath)
                    end_signal = Opti.algorithm_stats(Opti.population)
                    cwd = os.getcwd()
                    if Opti.postprocessing:
                        self.logger.info('Running Post-processing')
                        path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,0))
                        os.chdir(path)
                        if Opti.genealogytree:
                            pp.read_output(os.getcwd(),genealogytree=True,natoms=Opti.natoms)
                        else:
                            pp.read_output(os.getcwd(),genealogytree=False,natoms=Opti.natoms)
                        os.chdir(cwd)
                    if Opti.lattice_concentration:
                        if Opti.structure=='Defect':
                            self.logger.info('Running lattice concentration check')
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,rank))
                            os.chdir(path)
                            if Opti.BestIndsList:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+Opti.filename+'.xyz'))
                            else:
                                pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                            os.chdir(cwd)
                    Opti.write()
                    return True
        if Opti.generation !=0:
            self.logger.info("Mutating and evolving structures for generation %i" % Opti.generation)
        else:
            self.logger.info('Creating files for initial population')
        self.mutate_evolve_write_new_structures(Opti)
        self.logger.info("Changing status back.")
        self.flip_status_back()
        os.chdir(ingpath)
        Opti.write()
        return False
    
    def is_complete_looped(self):
        """Check if StructOpt looped generations are complete.
        """
        if 'filename' in self.structopt_parameters:
            opath = os.path.join(self.keywords['name'],self.structopt_parameters['filename'])
        else:
            opath = os.path.join(self.keywords['name'],'Output-rank0.txt')
        reachgrep=subprocess.Popen('grep "End of Execution" %s' % opath, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        reachrpt=reachgrep.communicate()[0]
        reachgrep.wait()
        if reachrpt=='':
            reachedend=False
        else:
            self.logger.info("StructOpt output file at %s shows calculation has finished." % opath)
            return True
        if 'filename' in self.structopt_parameters:
            optifile = os.path.join(os.path.join(self.keywords['name'], '{0}-rank0'.format(self.structopt_parameters['filename']),'Optimizer-restart-file.txt'))
        else:
            optifile = os.path.join(os.path.join(self.keywords['name'], 'Output-rank0','Optimizer-restart-file.txt'))
        ingpath = copy.deepcopy(self.keywords['name'])
        count = 0
        pwd = os.getcwd()
        while count < 10:
            if not os.path.isfile(optifile):
                os.chdir(ingpath)
                inputfile = os.path.join(ingpath,"structoptinput.txt")
                Opti = Optimizer(inputfile)
                Opti.algorithm_initialize()
                self.logger.info("Opti initialized")
                self.logger.info('Creating files for initial population')
                completeflag=True
            else:
                if not self.subfolders_complete_looped():
                    self.logger.info("Subfolders not complete. No further checks.")
                    self.flip_status_back()
                    completeflag = False
                else: 
                    self.logger.info("Evaluating completed structures.")
                    Opti = Optimizer(input=optifile)
                    self.logger.info('Optimizer read optimizer restart file')
                    self.logger.info("Generation %i" % Opti.generation)
                    Opti = self.evaluate_structures_for_fitness(Opti)
                    completeflag=True
                    if Opti.convergence:
                        os.chdir(ingpath)
                        end_signal = Opti.algorithm_stats(Opti.population)
                        if Opti.postprocessing:
                            self.logger.info('Running Post-processing')
                            path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,0))
                            os.chdir(path)
                            if Opti.genealogytree:
                                pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms)
                            else:
                                pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms)
                            os.chdir(cwd)
                        if Opti.lattice_concentration:
                            if Opti.structure=='Defect':
                                self.logger.info('Running lattice concentration check')
                                path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,rank))
                                os.chdir(path)
                                if Opti.BestIndsList:
                                    pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+Opti.filename+'.xyz'))
                                else:
                                    pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                                os.chdir(cwd)
                        Opti.write()
                        return True
            if completeflag:
                if Opti.generation !=0:
                    self.logger.info("Mutating and evolving structures for generation %i" % Opti.generation)
                else:
                    self.logger.info('Creating files for initial population')
                self.mutate_evolve_write_new_structures(Opti)
                self.logger.info("Changing status back.")
                self.flip_status_back()
                os.chdir(ingpath)
                Opti.write()
            print "Submitting from submission list."
            self.submit_from_submission_list()
            print "Clearing submission list."
            queue_commands.clear_submission_list()
#             mycwd=os.getcwd()
#             os.chdir(dirutil.get_mast_control_path())
#             mycommand=queue_commands.queue_submission_command("mastmon_submit.sh") #run the mastmon_submit.sh script in $MAST_CONTROL which creates a mastmon in order to do status checking on a compute node
#             mysub = subprocess.Popen(mycommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             mysub.wait()
#             os.chdir(mycwd)
            os.chdir(pwd)
            self.keywords['name'] = ingpath
            count +=1
            time.sleep(60)
        return False
    
    def is_ready_to_run(self):
        """Check if StructOpt ingredient is 
            ready to run.
        """
        dirname = self.keywords['name']
        flist = os.listdir(dirname)
        notready = 0
        if "structoptinput.txt" not in flist:
            notready += 1
        if self.potential_file:
            potfile = os.path.basename(self.potential_file)
            if potfile not in flist:
                notready += 1
        if "submit.sh" not in flist:
            notready = notready + 1
        if notready > 0:
            return False
        else:
            import sys
            ingpath = self.keywords['name']
            inputfile = os.path.join(ingpath,"structoptinput.txt")
            exec_line = self.keywords['program_keys']['mast_exec']
            sp_exec_line = exec_line.split(' ')
            new_exec_line = str()
            for segment in sp_exec_line:
                if ('.py' not in segment):
                    new_exec_line += segment + ' '
            try:
                optimizerpath = os.path.join(os.environ['MAST_CONTROL'],'Optimizer.py')
            except:
                raise MASTError(self.__class__.__name__,
                    "Cannot find Optimizer run file in 'structopt'! Current path: {0}".format(os.path.join(os.environ['MAST_CONTROL'],'Optimizer.py')))
            new_exec_line += ' {0} {1}'.format(optimizerpath,inputfile)
            self.keywords['program_keys']['mast_exec'] = new_exec_line
            return True
    
    def is_ready_to_next(self):
        if not self.subfolders_complete():
            self.logger.info("Subfolders not complete. No further checks.")
            self.flip_status_back()
            return False
        else:
            return True
    
    def is_started(self):
        """See if the ingredient has been started on
            the queue.
        """
        if os.path.isfile(os.path.join(self.keywords['name'],'structoptinput.txt')):
            return True
        else:
            return False
    
    def set_up_program_input(self):
        input_dict = self.structopt_parameters
        if ('restart_optimizer' in input_dict) and (input_dict['restart_optimizer']==True):
            if ('structopt_input_file' in input_dict) and \
                ('Optimizer-restart-file' in input_dict['structopt_input_file']):
                input_dict = self.transfer_restart_files(input_dict)
        if ('restart' in input_dict) and (input_dict['restart']==True):
            input_dict = self.transfer_restart_files(input_dict)
        self.structopt_parameters = input_dict
        if "MAST" in self.structopt_parameters['calc_method']:
            if 'LOOPED' in self.structopt_parameters['calc_method']:
                return self.set_up_program_input_looped()
            else:
                return self.set_up_program_input_one_gen()
        else:
            return self.set_up_program_input_whole()
    
    def set_up_program_input_whole(self):
        """Set up the StructOpt input file."""
        self.logger.info('Parameters for structopt: {0}'.format(self.structopt_parameters))
        if ('restart_optimizer' in self.structopt_parameters and self.structopt_parameters['restart_optimizer']):
            population_list = self.structopt_parameters['population']
            bests_list = self.structopt_parameters['BESTS']
        Optinit = Optimizer(self.structopt_parameters)
        ingpath = self.keywords['name']
        inputfile = os.path.join(ingpath,"structoptinput.txt")
        if 'potential_file' in self.structopt_parameters:
            newlocation = os.path.join(ingpath,os.path.basename(self.structopt_parameters['potential_file']))
            shutil.copyfile(self.structopt_parameters['potential_file'],newlocation)
            if not Optinit.pot_file:
                Optinit.pot_file = os.path.basename(self.structopt_parameters['potential_file'])
        if Optinit.restart_optimizer:
            Optinit.population = population_list
            Optinit.BESTS = bests_list
            Optinit.write(inputfile)
        else:
            Optinit.write(inputfile,restart=False)
        import sys
        exec_line = self.keywords['program_keys']['mast_exec']
        sp_exec_line = exec_line.split(' ')
        new_exec_line = str()
        for segment in sp_exec_line:
            if ('.py' not in segment):
                new_exec_line += segment + ' '
        try:
            structoptpath = os.path.join(os.environ['MAST_CONTROL'],'Optimizer.py')
        except:
            raise MASTError(self.__class__.__name__,
                "Cannot find Optimizer run file in 'structopt'! Current path: {0}".format(os.path.join(os.environ['MAST_CONTROL'],'Optimizer.py')))
        new_exec_line += ' {0} {1}'.format(structoptpath,inputfile)
        self.keywords['program_keys']['mast_exec'] = new_exec_line
        return
    
    def set_up_program_input_one_gen(self):
        """Set up the StructOpt input file."""
        Optinit = Optimizer(self.structopt_parameters)
        ingpath = self.keywords['name']
        inputfile = os.path.join(ingpath,"structoptinput.txt")
        if 'potential_file' in self.structopt_parameters:
            newlocation = os.path.join(ingpath,os.path.basename(self.structopt_parameters['potential_file']))
            shutil.copyfile(self.structopt_parameters['potential_file'],newlocation)
            if not Optinit.pot_file:
                Optinit.pot_file = os.path.basename(self.structopt_parameters['potential_file'])
        Optinit.write(inputfile,restart=False)
        self.keywords['program_keys']['mast_exec'] = "mast"
        return
    
    def set_up_program_input_looped(self):
        """Set up the StructOpt input file."""
        Optinit = Optimizer(self.structopt_parameters)
        ingpath = self.keywords['name']
        inputfile = os.path.join(ingpath,"structoptinput.txt")
        if 'potential_file' in self.structopt_parameters:
            newlocation = os.path.join(ingpath,os.path.basename(self.structopt_parameters['potential_file']))
            shutil.copyfile(self.structopt_parameters['potential_file'],newlocation)
            if not Optinit.pot_file:
                Optinit.pot_file = os.path.basename(self.structopt_parameters['potential_file'])
        Optinit.write(inputfile,restart=False)
        runmastpath = os.path.join(dirutil.get_mast_control_path(),"runmast.py")
        #mycommand=queue_commands.queue_submission_command("mastmon_submit.sh") 
        #"./"+os.path.join(mast_control,"mastmon_submit.sh")
        self.keywords['program_keys']['mast_exec'] = "python {0} >> $MAST_CONTROL/mastoutput 2> $MAST_CONTROL/errormast".format(runmastpath)
        return
      
    def eval_generation(self,optifile):
        self.logger.info("Evaluating completed structures.")
        Opti = Optimizer(input=optifile)
        self.logger.info('Optimizer read optimizer restart file')
        self.logger.info("Generation %i" % Opti.generation)
        Opti = self.evaluate_structures_for_fitness(Opti)
        #Check for convergence and run stats and postprocessing as desired if True
        if Opti.convergence:
            end_signal = Opti.algorithm_stats(Opti.population)
            if Opti.postprocessing:
                self.logger.info('Running Post-processing')
                path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,0))
                os.chdir(path)
                if Opti.genealogytree:
                    pp.read_output(os.getcwd(),genealogytree=True,natoms=self.natoms)
                else:
                    pp.read_output(os.getcwd(),genealogytree=False,natoms=self.natoms)
                os.chdir(cwd)
            if Opti.lattice_concentration:
                if Opti.structure=='Defect':
                    self.logger.info('Running lattice concentration check')
                    path = os.path.join(os.getcwd(), '{0}-rank{1}'.format(Opti.filename,rank))
                    os.chdir(path)
                    if Opti.BestIndsList:
                        pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'Bests-'+Opti.filename+'.xyz'))
                    else:
                        pp.get_lattice_concentration(os.path.join(os.getcwd(),'Bulkfile.xyz'),os.path.join(os.getcwd(),'indiv00.xyz'))
                    os.chdir(cwd)
            Opti.write()
            return
        self.logger.info("Mutating and evolving structures for generation %i" % Opti.generation)
        self.mutate_evolve_write_new_structures(Opti)
        self.logger.info("Changing status back.")
        self.flip_status_back()
        Opti.write()
        return
    
    def evaluate_structures_for_fitness(self, MyOpti):
        """Evaluate completed VASP runs' structures for fitness.
            Args:
                MyOpti <GA Optimizer>
            Returns:
                MyOpti <GA Optimizer>: optimizer with any 
                        settings that had been set
        """
        #Subfolders are complete. Now need to
        #evaluate them for fitness and restart the
        #Optimizer loop if necessary.
        population, offspring = self.get_complete_population(MyOpti)
        self.logger.info(os.environ)
        #evaluate fitness of new individuals in offspring
        for ind in offspring:
            outs = fitness_switch([MyOpti,ind])
            MyOpti.output.write(outs[1])
            ind=outs[0]
        for ind in offspring:
            self.logger.info("Offspring fitness: %3.3f" % ind.fitness)
            #make sure it worked.
        #Now append offspring to population, which includes
        #parents and fitnesses already.
        population.extend(offspring)
        self.logger.info("TTM Population: %s" % population)
        self.logger.info("Opti files: %s" % MyOpti.files)
        pwd = os.getcwd()
        ingpath = self.keywords['name']
        os.chdir(ingpath)
        population = MyOpti.generation_eval(population)
        os.chdir(pwd)
        self.logger.info("Generation incremented to %i" % MyOpti.generation)
        if MyOpti.convergence:
            self.logger.info("Converged!")
            #convokay = open("CONVERGED","wb")
            #convokay.write("Opti.convergence is True at %s\n" % time.asctime())
            #convokay.close()
            #end_signal = MyOpti.algorithm_stats(pop)
            return MyOpti
        else:
            self.logger.info("Clearing folders for next setup.")
            self.clear_folders()
        return MyOpti
    
    def get_complete_population(self, MyOpti):
        self.logger.info("Extracting VASP/LAMMPS outputs.")
        ingpath = self.keywords['name']
        population = list()
        offspring = list()
        pathtooutput = os.path.join(ingpath,MyOpti.filename+'-rank0')
        Totaloutputs = self.extract_output(ingpath)
        tolen = len(Totaloutputs)
        self.logger.info('Length of totaloutputs = {0}'.format(tolen))
        count = 0
        for idx in range(len(MyOpti.population)):
            ind = MyOpti.population[idx]
            if ind.energy != 0:
                population.append(ind)
                count = idx
            else:
                self.logger.info('Count = {0}'.format(count))
                try:
                    (myatoms, myenergy, mypressure) = Totaloutputs[idx-count]
                    passflag = True
                    self.logger.info('TTM atoms, energy, pressure: %s %s %s' % (myatoms, myenergy, mypressure))
                except:
                    passflag = False
                self.logger.info("TTM passflag: %s" % passflag)
                if passflag:
                    if MyOpti.structure=='Defect':
                        from MAST.structopt.tools import find_defects
                        outt=find_defects(myatoms.copy(), MyOpti.solidbulk, MyOpti.sf,
                            atomlistcheck=MyOpti.atomlist,trackvacs=MyOpti.trackvacs,
                            trackswaps=MyOpti.trackswaps,debug=False)
                        ind[0]=outt[0]
                        ind.buki=outt[1]
                        ind.vacancies = outt[2]
                        ind.swaps = outt[3]
                        MyOpti.output.write(outt[4])
                    elif MyOpti.structure=='Cluster':
                        myatoms.translate([-MyOpti.large_box_size/2.0,-MyOpti.large_box_size/2.0,-MyOpti.large_box_size/2.0])
                        ind[0] = myatoms.copy()
                    else:
                        ind[0] = myatoms.copy()
                    ind.energy = myenergy
                    ind.pressure = mypressure
                    offspring.append(ind)
        return population, offspring
    
    def extract_output(self, childname=""):
        """Extract output from the VASP/LAMMPS runs.
        """
        outputlist=list()
        keywords = copy.deepcopy(self.keywords)
        keywords['program_keys']['mast_exec'] = self.structopt_parameters['exec_mast']
        for subfolder in self.get_subfolder_list():
            keywords['name'] = os.path.join(childname,subfolder)
            if 'VASP' in self.structopt_parameters['calc_method']:
                mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
                mystrfile = os.path.join(subfolder,"CONTCAR")
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                mychecker = LammpsChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
                mystrfile = os.path.join(subfolder,"TRAJECTORY")
            mystructure = mychecker.get_structure_from_file(mystrfile)
            myatoms = AseAtomsAdaptor.get_atoms(mystructure)
            myenergy = mychecker.get_energy_from_energy_file()
            mypressure = mychecker.get_final_pressure()
            outputlist.append((myatoms, myenergy, mypressure))
        return outputlist

    def clear_folders(self):
        """Clear the Genetic Algorithm Evaluator ingredient.
            Assumes that the VASP ingredient has already
            been evaluated and the result passed
            to the child ingredient.
            Args:
                ingred_name <str>: Ingredient name to clear.
                    Match this carefully with a name from
                    the recipe.
        """
        fullpath = self.keywords['name']
        self.logger.info("Removing directories and files from ingredient specified at %s" % fullpath)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        os.mkdir(timestamp)
        for subfolder in self.get_subfolder_list():
            pth = os.path.join(fullpath, timestamp, os.path.basename(subfolder))
            os.renames(subfolder, pth)
        return
    
    def flip_status_back(self):
        """Flip status back to waiting."""
        ingpath = self.keywords['name']
        from MAST.ingredients.chopingredient import ChopIngredient
        cleared_ing = ChopIngredient(name=ingpath, program='None',program_keys = self.keywords['program_keys'], structure=self.keywords['structure'])
        cleared_ing.change_my_status("W")
        return

    def mutate_evolve_write_new_structures(self, MyOpti):
        """Mutate, evolve, and write new structures for VASP.
            Args:
                MyOpti <GA Optimizer>
            Returns:
                MyOpti <GA Optimizer>: optimizer with any 
                        settings that had been set
        """
        #Mutating and evolving
        ingpath = self.keywords['name']
        pathtooutput = os.path.join(ingpath,MyOpti.filename+'-rank0')
        offspring = MyOpti.generation_set(MyOpti.population)
        self.logger.info("generation set; pop len %i" % len(MyOpti.population))
        # Identify the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if ind.fitness==0]
        while len(invalid_ind) == 0:
            offspring = MyOpti.generation_set(offspring)
            invalid_ind = [ind for ind in offspring if ind.fitness==0]
        self.run_structures(invalid_ind,MyOpti)
        return MyOpti
    
    def run_structures(self, invalid_ind, MyOpti):
        ingpath = self.keywords['name']
        pathtooutput = os.path.join(ingpath,MyOpti.filename+'-rank0')
        #Delete all old files:
        dirlist = os.listdir(pathtooutput)
        for diritem in dirlist:
            if 'VASP' in self.structopt_parameters['calc_method']:
                if 'POSCAR_' in diritem:
                    os.remove(os.path.join(pathtooutput, diritem))
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                if 'DATA' in diritem:
                    os.remove(os.path.join(pathtooutput, diritem))
                if 'atom_symbols' in diritem:
                    os.remove(os.path.join(pathtooutput, diritem))
        #Write structures to POSCAR files or DATA files for evaluation in MAST
        MyOpti.output.write('length of individauls = {0}'.format(len(invalid_ind)))
        for inv_idx in range(len(invalid_ind)):
            MyOpti.output.write('This one = {0}'.format(inv_idx))
            if 'VASP' in self.structopt_parameters['calc_method']:
                indname = os.path.join(pathtooutput,'POSCAR_%02d' % inv_idx)
                #"%s/POSCAR_%02d" % (MyOpti.filename+'-rank0', i)
                import ase
                if MyOpti.structure == 'Defect':
                    indatoms = invalid_ind[inv_idx][0].copy()
                    indatoms.extend(invalid_ind[inv_idx].bulki)
                    #Run small optimizer to ensure atoms are not too close
                    #from ase.calculators.lj import LennardJones
                    #from ase.optimize import BFGS
                    #indatoms.set_calculator(LennardJones())
                    #dyn = BFGS(indatoms)
                    #dyn.run(fmax=0.01, steps=500)
                    #from MAST.structopt.tools.eval_energy import check_min_dist
                    min_len = 1.8 #Minimum distance between 2 atoms in angstroms
                    from ase.calculators.neighborlist import NeighborList
                    from ase import Atom, Atoms
                    import numpy
                    cutoffs=[2.0 for one in indatoms]
                    nl=NeighborList(cutoffs,bothways=True,self_interaction=False)
                    count = 0
                    dflag = True
                    while dflag:
                        count +=1
                        dflag = False
                        for one in indatoms:
                            nl.update(indatoms)
                            nbatoms=Atoms()
                            nbatoms.append(one)
                            indices, offsets=nl.get_neighbors(one.index)
                            for index, d in zip(indices,offsets):
                                index = int(index)
                                sym=indatoms[index].symbol
                                pos=indatoms[index].position + numpy.dot(d,indatoms.get_cell())
                                at=Atom(symbol=sym,position=pos)
                                nbatoms.append(at)
                            for i in range(1,len(nbatoms)):
                                d = nbatoms.get_distance(0,i)
                                if d < min_len:
                                    nbatoms.set_distance(0,i,min_len+.01,fix=0.5)
                                    dflag = True
                            idx = 0
                            for index, d in zip(indices,offsets):
                                index = int(index)
                                pos = nbatoms[idx+1].position - numpy.dot(d,indatoms.get_cell())
                                indatoms[index].position = pos
                                idx+=1
                            indatoms[one.index].position=nbatoms[0].position
                        if count > 1000:
                            self.logger.info("unable to satisfy min_len condition")
                            break
                    #indatoms, STR = check_min_dist(indatoms, MyOpti.structure, len(indatoms), min_len, '')
                    ase.io.write(indname, indatoms, "vasp", direct=True, sort=True, vasp5=True)
		    self.logger.info("HKK:: Writing individual {0}", indname) 
                else:
                    ase.io.write(indname,invalid_ind[i][0],"vasp", direct=True, sort=True, vasp5=True)
		    self.logger.info("HKK:: Writing invalid individual {0}", indname)
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                #indname = "%s/DATA_%02d" % (MyOpti.filename+'-rank0', i)
                indname = os.path.join(pathtooutput,'DATA_%02d' % i)
                if MyOpti.structure == 'Defect':
                    indatoms = invalid_ind[i][0]
                    indatoms.extend(invalid_ind[i].bulki)
                    inp_out.write_lammps_data(indname, indatoms)
                elif MyOpti.structure == 'Cluster':
                    indatoms = invalid_ind[i][0]
                    indatoms.set_cell([MyOpti.large_box_size,MyOpti.large_box_size,MyOpti.large_box_size])
                    indatoms.translate([MyOpti.large_box_size/2.0,MyOpti.large_box_size/2.0,MyOpti.large_box_size/2.0])
                    inp_out.write_lammps_data(indname, indatoms)
                else:
                    inp_out.write_lammps_data(indname, invalid_ind[i][0])
                shutil.copy(os.path.join(pathtooutput,'atom_symbols'), 
                    os.path.join(pathtooutput,'atom_symbols_{0:02d}'.format(i)))
#                 shutil.copy("{0}-rank0/atom_symbols".format(MyOpti.filename), 
#                     "{0}-rank0/atom_symbols_{1:02d}".format(MyOpti.filename,i))
            MyOpti.output.write(indname+'\n')
            #Append individuals in invalid_ind to population for writing out
            if MyOpti.generation > 0:
                MyOpti.population.append(invalid_ind[inv_idx])
        #Submit list of files to Tam's script to make subfolders
        self.make_subfolders_from_structures(pathtooutput)
        self.set_up_run_folders(pathtooutput)  #Start running.
        return
    
    def make_subfolders_from_structures(self, finddir="", childdir=""):
        """This method makes POSCAR/DATA-containing subfolders from POSCAR/DATA
            files found in a directory
            in childdir if a child directory is given, or in 
            mydir otherwise.
            Args:
                #mydir <str>: Ingredient directory (from self)
                finddir <str>: Directory containing list of structures
                childdir <str>: Child directory (optional)
            Files should be named:
                POSCAR_00 or DATA_00
                POSCAR_01 or DATA_01
                POSCAR_02 or DATA_02
                etc.
        """
        mydir = self.keywords['name']
        logger = logging.getLogger(mydir)
        logger = loggerutils.add_handler_for_recipe(mydir, logger)
        if childdir == "":
            childdir = mydir
        strfiles = os.listdir(finddir)
        for strfile in strfiles:
            if 'VASP' in self.structopt_parameters['calc_method']:
                if (not strfile[0:6] == "POSCAR"):
                    continue
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                if (not strfile[0:4] == "DATA"):
                    continue
            if "ase" in strfile:
                continue
            str_path = os.path.join(finddir, strfile)
            subname = strfile.split("_")[1]
            subpath = os.path.join(childdir, subname)
            if not os.path.isdir(subpath):
                os.mkdir(subpath)
            if 'VASP' in self.structopt_parameters['calc_method']:
                shutil.copy(str_path, "%s/POSCAR" % subpath)
                #Write Metadata
                from MAST.utility import Metadata
                metafile = Metadata(metafile='%s/metadata.txt' % subpath)
                metafile.write_data('directory created', time.asctime())
                metafile.write_data('name', subname)
                metafile.write_data('program', 'vasp')
                metafile.write_data('ingredient type', "BaseIngredient")
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                shutil.copy(str_path, "%s/DATA" % subpath)
                shutil.copy(os.path.join(finddir,'atom_symbols_{0}'.format(subname)),
                    "%s/atom_symbols" % subpath)
        return "Wrote all POSCAR/DATA files to subdirectories in %s" % childdir

    def set_up_run_folders(self, childname=""):
        """Set up the Optimizer ingredient subfolders
        """
        self.logger.info("setting up run folders")
        from MAST.ingredients.chopingredient import ChopIngredient
        keywords = self.keywords
        keywords['program_keys']['mast_exec'] = self.structopt_parameters['exec_mast']
        for subfolder in self.get_subfolder_list():
            keywords['name'] = os.path.join(childname,subfolder)
            if 'VASP' in self.structopt_parameters['calc_method']:
                keywords['program'] = 'vasp'
                mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                keywords['program'] = 'lammps'
                mychecker = LammpsChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            mychecker.set_up_program_input()
            mychoping = ChopIngredient(name=subfolder, program=keywords['program'], program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
            mychoping.write_submit_script()
            mychoping.run_singlerun()

    def get_subfolder_list(self):
        """Get subfolder list (numeric digits)
        """
        dircontents = os.listdir(self.keywords['name'])
        subfolders = list()
        for diritem in dircontents:
            fulldir = os.path.join(self.keywords['name'],diritem)
            if os.path.isdir(fulldir) and diritem.isdigit():
                subfolders.append(fulldir)
        subfolders.sort()
        return subfolders


    def subfolders_complete(self, childname=""):
        """Check if the Optimizer ingredient subfolders are complete
        """
        from MAST.ingredients.chopingredient import ChopIngredient
        subfolders = self.get_subfolder_list()
        if 'LAMMPS' in self.structopt_parameters['calc_method']:
            count = 2
        else:
            count = 1
        for check in range(count):
            allcomplete=0
            #HKK
            for subfolder in self.get_subfolder_list():
                if 'VASP' in self.structopt_parameters['calc_method']:
                    keywords = self.keywords
                    keywords['program'] = 'vasp'
                    mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
                elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                    keywords = self.keywords
                    keywords['program'] = 'lammps'
                    mychecker = LammpsChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
                if mychecker.is_complete():
                    allcomplete = allcomplete + 1
                    #print 'HKK :: Checking Subfolders,',subfolder, 'Complete'
                
                else:
          # HKK :: 04-20-15  :: CONTCAR was not updated when VASP run is incomplete. Fixed.
                    #print 'HKK :: Checking Subfolders,',subfolder, 'Incomplete'
                    #mychoping = ChopIngredient(name=subfolder, program=keywords['program'], program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
                    #mychoping.copy_file(copyfrom="CONTCAR", copyto="POSCAR") 
                    #mychoping.run_singlerun()
          # HKK :: 04-22-15  Stopping E evaluation if it iterated more than 3 times. Not going to converge. 
                    subdirlist = os.listdir(mychecker.keywords['name'])
                    count_opt_num=0
                    for file in subdirlist:
                        if any(i.isdigit() for i in file) is True:
                           count_opt_num = count_opt_num+1
                    path_oszicar = os.path.abspath('%s/VASP_replace/OSZICAR' % os.getenv("HOME"))
                    path_outcar = os.path.abspath('%s/VASP_replace/OUTCAR' % os.getenv("HOME"))
                    if count_opt_num >= 2:
                        print 'HKK :: Evaluated following structure more than 3 times. Copying low fitness files'
                        print 'HKK :: Checking Subfolders,',subfolder, 'Forced to complete'
                        shutil.copy(path_oszicar, subfolder)
                        shutil.copy(path_outcar, subfolder)
                        allcomplete = allcomplete + 1
                    else:
                        mychoping = ChopIngredient(name=subfolder, program=keywords['program'], program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
                        mychoping.copy_file_no_name_validation(copyfrom="CONTCAR", copyto="POSCAR")
                        print 'HKK :: Checking Subfolders,',subfolder, 'Incomplete'
                        mychoping.run_singlerun()
            if allcomplete == len(subfolders) and (allcomplete > 0):
                return True
            if 'LAMMPS' in self.structopt_parameters['calc_method']:
                time.sleep(60)
        return False
    
    def subfolders_complete_looped(self, childname=""):
        """Check if the Optimizer ingredient subfolders are complete
        """
        from MAST.ingredients.chopingredient import ChopIngredient
        subfolders = self.get_subfolder_list()
        allcomplete=0
        for subfolder in self.get_subfolder_list():
            if 'VASP' in self.structopt_parameters['calc_method']:
                keywords = self.keywords
                keywords['program'] = 'vasp'
                mychecker = VaspChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            elif 'LAMMPS' in self.structopt_parameters['calc_method']:
                keywords = self.keywords
                keywords['program'] = 'lammps'
                mychecker = LammpsChecker(name=subfolder, program_keys=self.keywords['program_keys'], structure=self.keywords['structure'])
            if mychecker.is_complete():
                allcomplete = allcomplete + 1
            else:
                mychoping = ChopIngredient(name=subfolder, program=keywords['program'], program_keys = self.keywords['program_keys'],structure=self.keywords['structure'])
                mychoping.run_singlerun()
        if allcomplete == len(subfolders) and (allcomplete > 0):
            return True
        return False
    
    def convert_asecalc2checker(self):
        return
    
    def get_structures(self, namelist):
        """Get structures from list
            Args:
                #mydir <str>: Ingredient directory (from self)
                namelist <list: List of file names of structures
            Returns:
                return_structures <list of Structure objects>
        """
        import numpy as np
        import pymatgen
        raise NotImplementedError()
        mydir = self.keywords['name']
        logger = logging.getLogger(mydir)
        logger = loggerutils.add_handler_for_recipe(mydir, logger)
        return_structures = list()
        
        for strline in namelist:
            tryline = strline.strip()
            if os.path.isfile(tryline):
                str_path = tryline
            else:
                trypath = os.path.join(mydir, tryline)
                if os.path.isfile(trypath):
                    str_path = trypath
                else:
                    logger.error("No file found at %s" % tryline)
                    continue
            if 'POSCAR' in str_path:
                one_structure = pymatgen.io.smartio.read_structure(str_path)
            elif 'xyz' in str_path:
                myxyz = pymatgen.io.xyzio.XYZ.from_file(str_path)
                allarr = np.array(myxyz.molecule.cart_coords)
                #Get box size:
                maxa=max(allarr[:,0])
                maxb=max(allarr[:,1])
                maxc=max(allarr[:,2])
                boxtol=0.01
                one_structure = myxyz.molecule.get_boxed_structure(maxa+boxtol,maxb+boxtol,maxc+boxtol)
            return_structures.append(one_structure)
        return return_structures

    def forward_final_structure_file(self, childpath, newname="Output-rank0"):
        """Forward the final structure.
            For StructOpt, this is the last population.
            Args:
                childpath <str>: Path of child ingredient
                newname <str>: new name (default 'OUTPUT')
        """
        if 'filename' in self.structopt_parameters:
            newname = self.structopt_parameters['filename']+'-rank0'
        path = os.path.join(childpath,newname)
        os.mkdir(path)
        files = os.listdir(os.path.join(os.getcwd(),newname))
        indivfiles = [one for one in files if ('indiv' in one) and ('iso' not in one)]
        self.logger.info('Found {0} individuals. Passing them to child folder at {1}'.format(len(indivfiles),path))
        for one in indivfiles:
            self.copy_a_file(path, newname+'/'+one, one)
        return
    
    def submit_from_submission_list(self):
        """Submit all entries from the submission list at
            $MAST_CONTROL/submitlist
            Adds a job number to the top of the "jobids" file in each
            ingredient directory.
        """
        control=dirutil.get_mast_control_path()
        submitlist=os.path.join(control, "submitlist")
        if not os.path.isfile(submitlist):
            print "No submission list at %s" % submitlist
            return
        submitfile=MASTFile(submitlist)
        #print 'HKK:: submitfile',submitfile
        subentries=list(submitfile.data)
        subentries.sort()
        submitted=dict()
        subcommand = "qsub submit.sh" #queue_commands.queue_submission_command()
        for subentry in subentries: # each entry is a directory name
            if len(subentry) == 0:
                continue
            subentry = subentry.strip()
            if len(subentry) == 0:
                continue
            if not os.path.isdir(subentry):
                continue
            elif subentry in submitted.keys():
                continue
            lastjobid = queue_commands.get_last_jobid(subentry)
            jstatus="not yet determined"
            if not (lastjobid == None):
                jstatus = queue_commands.get_job_status_from_queue_snapshot(subentry, lastjobid)
                if jstatus.lower() in ['r','q','h','e']: #running, queued, held, error
                    submitted[subentry]="Already on queue with status %s" % jstatus
                    continue
            os.chdir(subentry)
            subme=subprocess.Popen(subcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #subme.wait()
            status=subme.communicate()[0]
            queue_commands.write_to_jobids_file(subentry, status)
            submitted[subentry]=status
        queue_commands.print_submitted_dict(submitted)

# if __name__ == "__main__":
#     import sys
#     inputfile = sys.argv[1]
#     Opti = Optimizer(inputfile)
#     Opti.algorithm_initialize()
#     Opti.write()
#     self.mutate_evolve_write_new_structures(Opti)
