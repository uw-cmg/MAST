############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from dagutil import *
from MAST.utility import MASTError
class JobEntry(object):
    """JobEntry is an entry in class JobTable. \n
        isready : method to check if this job is ready or not \n
        sid : session id \n
        jid : job id \n
        name : jobname \n
        status : job status \n
        parents : set of parent jobs id \n
        completeparents : set of complete parent jobs \n
        ingredient_obj : ingredient obj \n
    """
    def __init__(self, sid, jid, jobname=None, type=None, ingredient=None):
        raise MASTError(self.__class__.__name__,"Do not use!")
        self.sid = sid # session id
        self.jid = jid # job id which is unique in a session.
        self.name = None # job name must be unique so that each job is stored in different directory
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
            self.name = ingredient.keywords['name']

        self.ingredient_obj = ingredient #ingredient object

        
    def addparent(self, jid):
        self.parents.add(jid)

    def addchild(self,jid):
        self.children.add(jid)
        
    def completeparent(self, jid):
        self.completeparents.add(jid)

    def is_ready(self):
        """is_ready returns whether this is ready or not. """
        return len(self.parents) == len(self.completeparents)
    
    def __str__(self):
        """Let me make this prettier later. """
        strpar = set2str(self.parents)
        strdonepar = set2str(self.completeparents)
        return '\t%d\t%d\t%s\t%d\t%s\t%s\t%s' % (self.sid, self.jid, self.name,
                                         self.status, self.type, strpar, strdonepar)
    def getinfo(self):
        """ getinfo returns list of members of JobEntry.
            
            See also:
                getheader, getformat, __str__
        """
        strpar = set2str(self.parents)
        strdonepar = set2str(self.completeparents)
        return [self.sid, self.jid, self.name,
                                         self.status, self.type, strpar, strdonepar]
    @classmethod
    def getheader(cls):
        """ getheader returns list of name of each field. This is a class method.
            
            See Also:
                getinfo, getformat, __str__
        """
        return ['sid','jid','jobname','status','type','parents','completeparents']

    @classmethod
    def getformat(cls):
        """ getformat returns a format to print out JobEntry information. This is a class method.
            
            See Also:
                getinfo, getheader, __str__
        """
        return "{:>5}{:>6}{:>40}{:>5}{:>20}{:>30}{:>20}"
    
