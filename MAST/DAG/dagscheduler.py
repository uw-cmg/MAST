import datetime
now = datetime.datetime.utcnow()

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
    for j in range(len(lA)):
        if j < len(lA)-1:
            elements += '{}, '.format(lA[j])
    elements += '{}'.format(lA[j])
    if len(elements) > maxlen-2:
        elements = elements[0:maxlen-5]+"..."
    out="{"+elements+"}"
    return out

Jobstatus = enum('PreQ','InQ','Complete')
JOB = Jobstatus
    
class JobEntry(object):
    def __init__(self, sid, jid, jobname, indir, outdir, type):
        self.sid = sid # session id
        self.jid = jid # job id which is unique in a session.
        self.name = jobname # job name may or may not be unique.
        self.indir = indir # input directory
        self.outdir = outdir # output directory
        self.status = JOB.PreQ # enum or integer
        self.type = type # string
        self.parents = set()
        self.completeparents = set()
        
    def addparent(self, jid):
        self.parents.add(jid)
        
    def completeparent(self, jid):
        self.completeparents.add(jid)

    def isready(self):
        return len(self.parents) == len(self.completeparents)
    
    def __str__(self):
        '''Let me make this prettier later '''
        strpar = str(self.parents)
        strdonepar = str(self.completeparents)
        return '\t%d\t%d\t%s\t%d\t%s\t%s\t%s' % (self.sid, self.jid, self.name,
                                         self.status, self.type, strpar, strdonepar)
    def getinfo(self):
        return 
                                   
class SessionEntry(object):

    def __init__(self, sid, totaljobs):
        self.sid =  sid
        self.totaljobs = totaljobs
        self.preq_jobs = totaljobs
        self.inq_jobs = 0
        self.completejobs = 0
        
    def submitjobs(njobs):
        self.preq_jobs -= njobs
        self.inq_jobs += njobs

    def completejobs(njobs):
        self.inq_jobs -= njobs
        self.complete_jobs += njobs

    def addnewjobs(njobs):
        self.totaljobs += njobs
        self.preq_jobs += njobs
        
    def iscomplete():
        return self.completejobs == self.totaljobs
    
    def __str__(self):
        return '\t%d\t%d\t%d\t%d\t%d' % (self.sid, self.totaljobs, self.preq_jobs,
                                         self.inq_jobs, self.completejobs)
    def getinfo(self):
        return [self.sid, self.totaljobs, self.preq_jobs,
                self.inq_jobs, self.completejobs]

    def getheader(self):
        return ['sid','total','preq','inq','complete']
        
    
class JobTable(object):

    def __init__(self):
        self.jobs = {}
        
    def addjob(self, job):
        if job.jid in self.jobs:
            raise Exception('JOB ID (jid=%d) CONFLICTED' % job.jid)
        self.jobs[job.jid] = job
        
    def deljob(self, jid):
        if jid not in self.jobs:
            raise Exception("JOB ID (jid=%d) DOESN'T EXIST" % jid)
        self.jobs.pop(jid)
        
class SessionTable(object):
    def __init__(self):
        self.sessions = {}
        
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
        header = s.getheader()
        row_format ="{:>5}{:>6}{:>5}{:>5}{:>9}"
        out = row_format.format(*header)+"\n"
        for row in data:
            out += row_format.format(*row)+"\n"
        return out

        

class Job(object):
    '''If previous version of dag scheduler or graph
        are not used, I will remove dictionary.
        I don't think we need dictionary for this if database
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
    ''' DAGPARSER has no member. So you don't need to make instance of this class.
    '''
    def __init__(self):
        pass

    @staticmethod
    def parse(jobs_file):
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
   
    @staticmethod
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
    def __init__(self):
        self.dag_parser    = DAGParser()
        pass
    
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

def main(jobs_file):
    sch_obj    = DAGScheduler()
    tasks_list = sch_obj.schedule(jobs_file)
    for index, tasks in enumerate(tasks_list):
        print "%d ==> %s\n" % (index + 1, ", ".join([job.plan() for job in tasks]))

if __name__ == "__main__":
    main("jobs_input1.txt")

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
