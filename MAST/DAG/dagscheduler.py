"""
This is for schedulring.
@author: Hyunwoo Kim
"""
import datetime
import time

now = datetime.datetime.utcnow()
MAXSID = 10000
MAXJID = 100000
SCHEDULING_INTERVAL = 10

def enum(**enums):
    return type('Enum',(),enums)

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def set2str(A,maxlen=20):

    if type(A) is not set:
        print 'INPUT IS NOT A SET.'
        return ""
    
    out = "{"
    lA = list(A)
    elements = ""
    nelem = len(lA)
    for j in range(nelem):
        if j < len(lA)-1:
            elements += '{}, '.format(lA[j])

    if nelem >0:
        elements += '{}'.format(lA[nelem-1])
        
    if len(elements) > maxlen-2:
        elements = elements[0:maxlen-5]+"..."
    out="{"+elements+"}"
    return out

Jobstatus = enum('PreQ','InQ','Complete','Error')
JOB = Jobstatus
    
class JobEntry(object):
    """JobEntry is an entry in class JobTable.
        isready : method to check if this job is ready or not
        sid : session id
        jid : job id
        name : jobname
        status : job status
        parents : set of parent jobs id
        completeparents : set of complete parent jobs
        ingredient_obj : ingredient obj
    """
    def __init__(self, sid, jid, jobname=None, type=None, ingredient=None):
        self.sid = sid # session id
        self.jid = jid # job id which is unique in a session.
        self.name = jobname # job name may or may not be unique.
        self.status = JOB.PreQ # enum or integer
        self.parents = set() # jid of parents
        self.completeparents = set() #jid of complete parents
        self.children = set() # jid of children
        
        if ingredient is not None:
            self.type = gettype(ingredient)
        else:
            self.type = type # string

        if jobname is not None:
            self.name = jobname
        elif ingredient is not None:
            self.name = ingredient.keywords['name'].split('/')[-1]
        else:
            self.name = 'noname'

        self.ingredient_obj = ingredient #ingredient object

        
    def addparent(self, jid):
        self.parents.add(jid)

    def addchild(self,jid):
        self.children.add(jid)
        
    def completeparent(self, jid):
        self.completeparents.add(jid)

    def is_ready(self):
        '''isready returns whether this is ready or not'''
        return len(self.parents) == len(self.completeparents)
    
    def __str__(self):
        '''Let me make this prettier later '''
        strpar = set2str(self.parents)
        strdonepar = set2str(self.completeparents)
        return '\t%d\t%d\t%s\t%d\t%s\t%s\t%s' % (self.sid, self.jid, self.name,
                                         self.status, self.type, strpar, strdonepar)
    def getinfo(self):
        strpar = set2str(self.parents)
        strdonepar = set2str(self.completeparents)
        return [self.sid, self.jid, self.name,
                                         self.status, self.type, strpar, strdonepar]
    @classmethod
    def getheader(cls):
        return ['sid','jid','jobname','status','type','parents','completeparents']

    @classmethod
    def getformat(cls):
        return "{:>5}{:>6}{:>8}{:>5}{:>5}{:>15}{:>20}"
    
class SessionEntry(object):

    def __init__(self, sid, totaljobs):
        self.sid =  sid
        self.totaljobs = totaljobs
        self.preq_jobs = totaljobs
        self.inq_jobs = 0
        self.completejobs = 0
        
    def set_jobs_to_submitted(self, njobs):
        self.preq_jobs -= njobs
        self.inq_jobs += njobs

    def set_jobs_to_complete(self, njobs):
        self.inq_jobs -= njobs
        self.completejobs += njobs

    def add_newjobs(self, njobs, status = JOB.PreQ):
        self.totaljobs += njobs
        if status is JOB.PreQ:
            self.preq_jobs += njobs
        elif status is JOB.InQ:
            self.inq_jobs += njobs
        elif status is JOB.Complete:
            self.completejobs += njobs
        
    def is_complete(self):
        return self.completejobs == self.totaljobs
    
    def __str__(self):
        return '\t%d\t%d\t%d\t%d\t%d' % (self.sid, self.totaljobs, self.preq_jobs,
                                         self.inq_jobs, self.completejobs)
    def getinfo(self):
        return [self.sid, self.totaljobs, self.preq_jobs,
                self.inq_jobs, self.completejobs]

    @classmethod
    def getheader(cls):
        return ['sid','total','preq','inq','complete']
    
    @classmethod
    def getformat(cls):
        return "{:>5}{:>6}{:>5}{:>5}{:>9}"
    
class JobTable(object):

    def __init__(self):
        self.jobs = {}
        self._nextjid = 1
        self._name2jid = {}

    def add_dependency(self, parent, child):
        self.jobs[self._name2jid[child]].addparent(self._name2jid[parent])
        self.jobs[self._name2jid[parent]].addchild(self._name2jid[child])

    def addjob(self, job):
        if job.jid in self.jobs:
            raise Exception('JOB ID (jid=%d) CONFLICTED' % job.jid)
        self.jobs[job.jid] = job
        self._name2jid[job.name] = job.jid

    def deljob(self, jid):
        if jid not in self.jobs:
            raise Exception("JOB ID (jid=%d) DOESN'T EXIST" % jid)
        self.jobs.pop(jid)

    def __str__(self):
        import numpy as np
        lists = []
        for j in self.jobs.itervalues():
            lists.append(j.getinfo())
        data = np.array(lists)

        header = JobEntry.getheader()
        row_format = JobEntry.getformat()
        out = row_format.format(*header)+"\n"
        for row in data:
            out += row_format.format(*row)+"\n"
        return out
    
    def __len__(self):
        '''len(jobtableobj) = number of jobs'''
        return len(self.jobs)

    def get_jid(self):
        '''get_sid returns unique session id in a session table object.'''
        while (self._nextjid in self.jobs):
            self._nextjid = (self._nextjid ) % MAXJID + 1
        return self._nextjid
    
    def update_complete_parent_set(self, children_jids, complete_parent_jid):
        if type(children_jids) is not list:
            children_jids = list(children_jids)
        for cjid in children_jids:
            self.jobs[cjid].completeparent(complete_parent_jid)
            
class SessionTable(object):
    def __init__(self):
        self.sessions = {}
        self._nextsid = 1
        
    def completejobs(self, sid, njobs):
        self.sessions[sid].set_jobs_to_complete(njobs)
        
    def submitjobs(self, sid, njobs):
        self.sessions[sid].set_jobs_to_submitted(njobs)

    def runjobs(self, sid, njobs):
        self.submitjobs(sid, njobs)
        
    def addsession(self, session):
        if session.sid in self.sessions:
            raise Exception('SESSION ID (sid=%d) CONFLICTED' % session.sid)
        self.sessions[session.sid] = session
        
    def delsession(self, sid):
        if sid not in self.session:
            raise Exception("SESSION ID (sid=%d) DOESN'T EXIST" % sid)
        self.sessions.pop(sid)

    def __str__(self):
        import numpy as np
        lists = []
        for s in self.sessions.itervalues():
            lists.append(s.getinfo())
        data = np.array(lists)
        header = SessionEntry.getheader()
        row_format = SessionEntry.getformat()
        out = row_format.format(*header)+"\n"
        for row in data:
            out += row_format.format(*row)+"\n"
        return out
    
    def get_sid(self):
        '''get_sid returns unique session id in a session table object.'''
        while (self._nextsid in self.sessions):
            self._nextsid = (self._nextsid ) % MAXSID + 1
        return self._nextsid

    def addnewjobs(self, sid, njobs, status = JOB.PreQ):
        self.sessions[sid].addnewjobs(njobs,status)
        
class Job(object):
    '''If previous version of dag scheduler or graph
        are not used, I will remove dictionary.
        I donot think we need dictionary for this if database
        is used for scheduling. hwkim
    '''
    def __init__(self, name):
        self.name = name
        self.jid = '' #Some jobs may have same name.
        self.pre_scripts = [] 
        self.post_scripts = [] #let's discuss about post and pre script
        self.parents_name = [] #potentially deprecated
        self.parents_set = set()
        
    def add_parents_name(self,pname):
        self.parents_name.append(pname)
        
    def add_parents_id(self,pid):
        self.parents_name.append(pid)
        
    def add_pre_script(self, script):
        self.pre_scripts.append(script)

    def add_post_script(self, script):
        self.post_scripts.append(script)

    def plan(self):
        plan = "( "
        if self.pre_scripts:
            plan += ", ".join(self.pre_scripts)
            plan += " -> "
        plan += self.name
        if self.post_scripts:
            plan += " -> "
            plan += ", ".join(self.post_scripts)
        plan += " )"
        return plan
    
    def get_jobinfo():
        info ={ "job_id":self.id,
                "job_name":self.name,
                "start": None,
                "end": None,
                "session_id": None
            }
        return info

class IDManager(object):
    def __init__(self, jobtable, sessiontable):
        self.nextjid = 0
        self.nexsid = 0
        self.jobtable = jobtable
        self.sessiontable = sessiontable
        
    def get_jid(self):
        self.nextjid = self.nextjid+1
        while self.nextjid in self.jobtable.jidset:
            self.nextjid = self.nextjid+1
        return self.nextjid
    
    def get_sid(self):
        self.nextsid = self.nextsid+1
        while self.nextsid in self.jobtable.sidset:
            self.nextsid = self.nextsid+1
        return self.nextsid
    
class DAGParser:
    ''' DAGPARSER has no member. So you do not need to make instance of this class.
    '''
    def __init__(self):
        pass

    @classmethod
    def parse(cls, jobs_file):
        '''Parses the jobs file and creates a dictionary, 
           with parent job as key and the dependent jobs as
           children.
           JOBS_FILE is input file name. This should be modified later for recipe template or
           something else.
        '''
        dependency_dict = {}
        jobs_dict       = {}
        f_ptr    = open(jobs_file, "r")

        for line in f_ptr.readlines():
            line = line.strip()
            #validate input line
            if not line or line.startswith('#'):
                continue

            #line elements
            elts = [elt.strip() for elt in line.split()]

            #Actions based on initial keyword
            init_keyword = elts[0].lower()

            #Job Keyword
            if init_keyword == "ingredient": 
                job_obj   = Job(elts[1])
                jobs_dict.setdefault(elts[1], job_obj)
                dependency_dict.setdefault(job_obj, set())
                continue

            #script keyword
            if init_keyword == "script":
                if elts[1].lower() == "pre":
                    job_obj   = jobs_dict[elts[2]] 
                    job_obj.add_pre_script(elts[3])
                elif elts[1].lower() == "post":
                    job_obj   = jobs_dict[elts[2]]
                    job_obj.add_post_script(elts[3])
                continue

            #parent keyword
            if init_keyword == "parent":
                parent_objs = []
                is_parent   = True
                for elt in elts[1:]:
                    if elt.lower() == "child":
                        is_parent = False
                        continue
                    
                    if is_parent:
                        parent_objs.append(jobs_dict[elt])

                    else:
                        for parent in parent_objs:
                            dependency_dict[jobs_dict[elt]].add(parent)
            	#parent child set
		(pset,cset) = get_parent_child_sets(elts)
		#add parent jobs to child jobs 
	
        f_ptr.close()
        return dependency_dict
   
    @classmethod
    def get_jobinfo(jobs_file):
	'''GET_JOBINFO reads input files and returns dictionary of jobs; jobs_dict. 
	    Each job has dependency, job name. This assumes that everything is a job. 
	'''
	dependency_dict = {}
        jobs_dict       = {}
        f_ptr    = open(jobs_file, "r")

        for line in f_ptr.readlines():
            line = line.strip()
            #validate input line
            if not line or line.startswith('#'):
                continue

            #line elements
            elts = [elt.strip() for elt in line.split()]

            #Actions based on initial keyword
            init_keyword = elts[0].lower()

            #Job Keyword
            if init_keyword == "ingredient": #TTM 2013-03-19 change init_keyword "job" to "ingredient"
                job_obj   = Job(elts[1])
                jobs_dict.setdefault(elts[1], job_obj)
                dependency_dict.setdefault(job_obj, set())
                continue

            #script keyword
            if init_keyword == "script":
                if elts[1].lower() == "pre":
                    job_obj   = jobs_dict[elts[2]] 
                    job_obj.add_pre_script(elts[3])
                elif elts[1].lower() == "post":
                    job_obj   = jobs_dict[elts[2]]
                    job_obj.add_post_script(elts[3])
                continue

            #parent keyword
            if init_keyword == "parent":
            	#parent child set
		(pset,cset) = get_parent_child_sets(elts)
		# safeguard for post and pre script
		pset = set.intersection(pset,set(jobs_dict.keys()))
		cset = set.intersection(cset,set(jobs_dict.keys()))
		#add parent jobs to child jobs 
	    for j in list(cset):
		jobs_dict[j].parents_set.update(pset)
	
        f_ptr.close()
        return jobs_dict 

class DAGScheduler:
    """In a session, the maximum number of jobs is MAXJID at once."""
    def __init__(self):
        self.jobtable = JobTable()
        self.sessiontable = SessionTable()
        self._run_mode = 'noqsub' # noqsub, serial (by qsub), parallel

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

    def addjobs(self, ingredients_dict, dependency_dict, sid=None):
        '''If sid is not specified, then new unique sid will be assigned.'''

        if sid is None:
            sid = self.sessiontable.get_sid()

        njobs = len(ingredients_dict)
        # Add a new session 
        if not self.has_session(sid):
            asession = SessionEntry(sid,njobs)
            self.sessiontable.addsession(asession)
        else:
            self.sessiontable[sid].addnewjobs(njobs)
            
        # Insert jobs
        for name, ingredient in ingredients_dict.iteritems():
            jid = self.jobtable.get_jid()
            ajob = JobEntry(sid=sid, jid=jid, jobname=name, ingredient=ingredient)
            print jid, name
            self.jobtable.addjob(ajob)
            
        # Update depdict
        for child, parent in dependency_dict.iteritems():
            for aparent in parent:
                # For DEBUGGING
                print child +"<="+ aparent
                self.jobtable.add_dependency(aparent,child)

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
            time.sleep(SCHEDULING_INTERVAL)
            
    def _run(self):
        print ':: run jobs ::'
        self.run_jobs()
        print ':: update job status ::'
        self.update_job_status()
        
        
# This part is for developing job and parser
def get_parent_child_sets(tokens):
    '''This function returns two sets as a tuple.
        Ex) (pset, cset) = get_parent_child_sets(tokens):
    '''
    parentset = set()
    childset = set()

    isparent = False
    ischild = False
    for item in tokens:
        if item.lower() == 'parent':
            isparent = True
            continue
        elif item.lower() == 'child':
            ischild = True
            isparent = False
            continue

        if isparent:
            parentset.add(item)
        if ischild:
            childset.add(item)
    return (parentset, childset)

def gettype(ingredientobj):
    stype=str(type(ingredientobj))
    stype = stype[1:-1]
    tokens = stype.replace("'","").split('.')
    return tokens[-1]
