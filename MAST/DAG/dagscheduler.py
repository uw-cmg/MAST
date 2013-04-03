import datetime
now = datetime.datetime.utcnow()

class Job(object):
    '''If previous version of dag scheduler or graph
        are not used, I will remove dictionary.
        I don't think we need dictionary for this if database
        is used for scheduling. hwkim
    '''
    def __init__(self, name):
        self.name = name
        self.id = '' #for speed ?
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
        return 

class Session(object):
    '''This is a set of job. Jobs will be scheduled in a session.
        All jobs out of a session are independent.
    '''
    def __init__(self,name,id):
        self.name = name
        self.id = id #unique
        self.jobs = {} #list of jobs
        self.start = None #datetime
        self.end = None #datetime
        
    def addjob(job):
        self.jobs.append(job)
        
    def getjobsid():
        jids = []
        for job in self.jobs:
            jids.append(job.id)
        return jids
    
    def getinfo():
        jids = self.jobs.getjobsid()
        info = {"session_id":"d_1",
                "jobs_id":jids,
                "start":now(),
                "end":None
            }
        return info

class DAGParser:
    ''' DAGPARSER has no member. So you don't need to make instance of this class. '''
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

