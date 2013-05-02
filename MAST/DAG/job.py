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
