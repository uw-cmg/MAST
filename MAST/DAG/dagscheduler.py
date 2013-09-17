from dagutil import *
from jobentry import JobEntry
from jobtable import JobTable
from sessionentry import SessionEntry
from sessiontable import SessionTable
import subprocess
import os
from MAST.utility import MASTError

"""
This is for schedulring.
@author: Hyunwoo Kim
"""

class DAGScheduler:
    """In a session, the maximum number of jobs is MAXJID at once."""
    
    def __init__(self):
        raise MASTError(self.__class__.__name__,"Do not use!")
        self.jobtable = JobTable()
        self.sessiontable = SessionTable()
        self._run_mode = 'noqsub' # noqsub, serial (by qsub), parallel
        self.home = os.path.expandvars(os.environ['MAST_SCRATCH'])
        self._ARCHIVE = os.path.expandvars(os.environ['MAST_ARCHIVE'])

    def set_mode(self,mode):
        self._run_mode = mode
        
    def has_complete_session(self):
        '''check all sessions in session table and returns number of complete session'''
        ncomplete_sessions=0
        for asession in self.sessiontable.sessions.itervalues():
            if asession.is_complete():
                ncomplete_sessions = ncomplete_sessions +1
        return ncomplete_sessions
    
    def get_complete_session_id(self):
        """returns complete session id (sid)"""
        csid = []
        for asession in self.sessiontable.sessions.itervalues():
            if asession.is_complete():
                csid.append(asession.sid)
        return csid
    
    def has_incomplete_session(self):
        return not len(self.sessiontable.sessions) == self.has_complete_session()
        
    def has_session(self, sid):
        return sid in self.sessiontable.sessions

    def addjobs(self, ingredients_dict, dependency_dict,sname=None, sid=None):
        '''If sid is not specified, then new unique sid will be assigned.'''

        if sid is None:
            sid = self.sessiontable.get_new_sid()

         # Add a new session 
        if not self.has_session(sid):
            asession = SessionEntry(sid, 0, sname) 
            self.sessiontable.addsession(asession)

        session_dir = os.path.join(self.home, sname)
        # Insert jobs
        for name, ingredient in ingredients_dict.iteritems():
            jid = self.jobtable.get_jid()
            self.sessiontable.add_newjob(sid, jid)

            # attach session_dir to name of ingredient
            ingredient.keywords['name'] = os.path.join(session_dir, ingredient.keywords['name'])
            # attach session_dir to name of child ingredients
            child_dict = ingredient.keywords['child_dict']
            keys = child_dict.keys()
            for oldkey in keys:
                newkey = os.path.join(session_dir, oldkey)
                child_dict[newkey] = child_dict[oldkey]
                del child_dict[oldkey]

            ajob = JobEntry(sid=sid, jid=jid, ingredient=ingredient)
            self.jobtable.addjob(ajob)
            
        # Update depdict
        for child, parent in dependency_dict.iteritems():
            for aparent in parent:
                # For DEBUGGING
                #print child +"<="+ aparent
                pname = os.path.join(session_dir,aparent)
                kidname = os.path.join(session_dir,child)
                #print kidname +"<="+ pname
                self.jobtable.add_dependency(pname,kidname)

    def update_job_status(self, verbose=0):
        '''udpate session table and update children'''
        for jid, job in self.jobtable.jobs.iteritems():
            if job.status == JOB.CompWait:
                if verbose == 1:
                    print "COMPWAIT job.name: ",job.name
                    print "    job.status: ",job.status
                    print "    job.is_ready():",job.is_ready()
                    print "    job.ingredient_obj.is_ready_to_run():",job.ingredient_obj.is_ready_to_run()
                self.jobtable.jobs[jid].status = JOB.Complete
                if verbose == 1:
                    print "    new job.status: ", self.jobtable.jobs[jid].status
                self.sessiontable.completejobs(job.sid, njobs=jid)
                job.ingredient_obj.update_children()
                self.jobtable.update_complete_parent_set(job.children,jid)
                
    def run_jobs(self, verbose=0):
        for jid, job in self.jobtable.jobs.iteritems():

            if verbose == 1:
                print "job.name: ",job.name
                print "    job.status: ",job.status
                print "    job.is_ready():",job.is_ready()
                print "    job.ingredient_obj.is_ready_to_run():",job.ingredient_obj.is_ready_to_run()
            
            if job.status == JOB.Complete:
                pass
            
            elif job.ingredient_obj.is_complete():
                if verbose == 1:
                    print "    job.ingredient_obj.is_complete(): True"
                job.status = JOB.Complete #TTM go straight to complete? Without running put jobs into InQ. In next turn, it will be got complete status
                
                self.sessiontable.runjobs(sid=job.sid, njobs=jid) # update session table
                if verbose == 1:
                    print "    new job.status: ", self.jobtable.jobs[jid].status
                self.sessiontable.completejobs(job.sid, njobs=jid)
                job.ingredient_obj.update_children()
                self.jobtable.update_complete_parent_set(job.children,jid)
                
            elif job.status is JOB.PreQ and job.is_ready() and job.ingredient_obj.is_ready_to_run():
                job.ingredient_obj.run() #TTM do not default mode in dag scheduler; make ingredients default. mode=self._run_mode)
                job.status = JOB.InQ
                sid = job.sid
                self.sessiontable.runjobs(sid=sid, njobs=jid) # update session table
                
            elif job.status is JOB.PreQ and job.is_ready():
                job.ingredient_obj.write_files()
            
            #elif self.jobtable.jobs[jid].status == JOB.CompWait:
            #    if verbose == 1:
            #        print "COMPWAIT job.name: ",job.name
            #        print "    job.status: ",job.status
            #        print "    job.is_ready():",job.is_ready()
            #        print "    job.ingredient_obj.is_ready_to_run():",job.ingredient_obj.is_ready_to_run()
            #    self.jobtable.jobs[jid].status = JOB.Complete
            #    if verbose == 1:
            #        print "    new job.status: ", self.jobtable.jobs[jid].status
            #    self.sessiontable.completejobs(job.sid, njobs=jid)
            #    job.ingredient_obj.update_children()
            #    self.jobtable.update_complete_parent_set(job.children,jid)

        
    def schedule(self, jobs_file):
        '''Parses the input file and applies topological sorting
           to find the set of jobs in topological order
        '''
        dependency_dict = self.dag_parser.parse(jobs_file)
        tasks_list      = self.sort_tasks(dependency_dict)
        return tasks_list

    def sort_tasks(self, dependency_dict):
        '''Apply topological sorting to fetch a list, which has
           sets of tasks that can be scheduled in parallel
        '''
        tasks_list = []
        while True:
            independent_jobs = set(parent for parent, dep_jobs in dependency_dict.iteritems() if not dep_jobs)
            if not independent_jobs:
                break
            tasks_list.append(independent_jobs)
            dependency_dict  = { parent : (dep_jobs - independent_jobs) for parent, dep_jobs in dependency_dict.iteritems() if parent not in independent_jobs}

        if dependency_dict:
             print "Cyclic dependencies exist between %s !!!" ", ".join(repr(entry) for entry in dependencies.iteritems())
             return []
        return tasks_list

    def show_session_table(self):
        print ' '
        print self.sessiontable

    def run(self, niter=None, verbose=0):
        """Run DAG scheduler.
            Args:
                niter <int>: number of iterations to run
                            None (default) = run until all sessions are complete
                            Otherwise, run scheduler niter number of times
                verbose <int>: 0 (defult) = not verbose
                               1 = include extra notification messages
            Returns:
                csnames <set of str>: set of completed session names; path of
                    complete sessions relative to MAST_SCRATCH directory
        """
        iter = 0
        csnames = set()
        while self.has_incomplete_session():
            self._run(verbose=verbose)
            #TTM DEBUG: do not move sessions to archive directory within loop
            iter = iter+1
            if niter is not None and iter >= niter:
                #print "TTM DEBUG: break"
                break
            time.sleep(SCHEDULING_INTERVAL) #TTM move sleep

        csnames = csnames.union(self._move_to_archive())
        return csnames

            
    def _run(self, verbose=0):
        self.run_jobs(verbose)
        #self.update_job_status(verbose)
        
    def _get_session_path(self, sid):
        return os.path.join(self.home,self.sessiontable.get_sname(sid))
    
    def del_session(self,sid):
        #print self.sessiontable
        self.sessiontable.delsession(sid)
        self.jobtable.del_session(sid)

    def _move_to_archive(self):
        csid = self.get_complete_session_id()
        os.chdir(self.home)
        dst = self._ARCHIVE
        csnames = set()
        if not os.path.exists(dst):
            os.mkdir(dst)
            
        for sid in csid:
            src = self._get_session_path(sid)
            os.system("mv %s %s"% (src,dst))
            csnames.add(self.sessiontable.get_sname(sid))
            self.del_session(sid)
            
        return csnames

    def __repr__(self):
        return str(self.sessiontable)            
