import datetime
now = datetime.datetime.utcnow()

class Job(object):
    def __init__(self, name):
        self.name = name
        self.id = ''
        self.pre_scripts = []
        self.post_scripts = []
        self.parents_name = []
        self.parents_id = [] #for unique jobid in DB

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
                "dish_id": None
            }
        return 

class Dish(object):
    def __init__(self,name,id):
        self.name = name
        self.id = id #unique
        self.jobs = [] #list of jobs
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
        info = {"dish_id":"d_1",
                "jobs_id":jids,
                "start":now(),
                "end":None
            }
        return info
    
class DAGParser:
    def __init__(self):
        pass

    def parse(self, jobs_file):
        '''Parses the jobs file and creates a dictionary, 
           with parent job as key and the dependent jobs as
           children
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
            if init_keyword == "job":
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

        f_ptr.close()
        return dependency_dict
                    

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
