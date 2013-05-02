from dagutil import *

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
    
