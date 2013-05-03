from dagutil import *
from jobentry import JobEntry
from jobtable import JobTable
from sessionentry import SessionEntry
from sessiontable import SessionTable
import subprocess
import os
"""
This is for schedulring.
@author: Hyunwoo Kim
"""

class DAGScheduler:
    """In a session, the maximum number of jobs is MAXJID at once."""
    def __init__(self):
        self.jobtable = JobTable()
        self.sessiontable = SessionTable()
        self._run_mode = 'noqsub' # noqsub, serial (by qsub), parallel
        self.home = os.environ['MAST_SCRATCH']

    def set_mode(self,mode):
        self._run_mode = mode
        
    def has_complete_session(self):
        '''check all sessions in session table and returns number of complete session'''
        ncomplete_sessions=0
        for assession in self.sessiontable.sessions.itervalues():
            if assession.is_complete():
                ncomplete_sessions = ncomplete_sessions +1
        return ncomplete_sessions
    
    def has_incomplete_session(self):
        return not len(self.sessiontable.sessions) == self.has_complete_session()
        
    def has_session(self, sid):
        return sid in self.sessiontable.sessions

    def addjobs(self, ingredients_dict, dependency_dict,sname=None, sid=None):
        '''If sid is not specified, then new unique sid will be assigned.'''

        if sid is None:
            sid = self.sessiontable.get_sid()

        njobs = len(ingredients_dict)
        # Add a new session 
        if not self.has_session(sid):
            asession = SessionEntry(sid, njobs, sname)
            self.sessiontable.addsession(asession)
        else:
            self.sessiontable[sid].addnewjobs(njobs)

        session_dir = os.path.join(self.home, sname)
        # Insert jobs
        for name, ingredient in ingredients_dict.iteritems():
            jid = self.jobtable.get_jid()
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
            print jid, name
            self.jobtable.addjob(ajob)
            
        # Update depdict
        for child, parent in dependency_dict.iteritems():
            for aparent in parent:
                # For DEBUGGING
                print child +"<="+ aparent
                pname = os.path.join(session_dir,aparent)
                kidname = os.path.join(session_dir,child)
                print kidname +"<="+ pname
                self.jobtable.add_dependency(pname,kidname)

    def update_job_status(self):
        '''udpate session table and update children'''
        for jid, job in self.jobtable.jobs.iteritems():
            if job.status == JOB.InQ and job.ingredient_obj.is_complete():
                job.status = JOB.Complete
                print '\njob [%s] is complete' % job.name
                self.sessiontable.completejobs(job.sid, njobs=1)
                job.ingredient_obj.update_children()
                self.jobtable.update_complete_parent_set(job.children,jid)
                
    def run_jobs(self):
        for jid, job in self.jobtable.jobs.iteritems():
            if job.status is JOB.PreQ and job.is_ready() and job.ingredient_obj.is_ready_to_run():
                print 'run [sid=%d,jid=%d] %s' % (job.sid, job.jid, job.name)
                job.ingredient_obj.run() #TTM do not default mode in dag scheduler; make ingredients default. mode=self._run_mode)
                job.status = JOB.InQ
                sid = job.sid
                self.sessiontable.runjobs(sid=sid, njobs=1) # update session table
                
            if job.status is JOB.PreQ and job.is_ready():
                print 'write files [sid=%d,jid=%d] %s' % (job.sid, job.jid, job.name)
                job.ingredient_obj.write_files()

        
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

    def run(self):
        while self.has_incomplete_session():
            self._run()
            print 'I am at %s' % os.getcwd()
            time.sleep(SCHEDULING_INTERVAL)
            
    def _run(self):
        print ':: run jobs ::'
        self.run_jobs()
        print ':: update job status ::'
        self.update_job_status()
        
        
