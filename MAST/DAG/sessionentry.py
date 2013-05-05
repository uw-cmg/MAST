from dagutil import *

class SessionEntry(object):

    def __init__(self, sid, totaljobs, name=None):
        self.sid =  sid
        self.totaljobs = totaljobs
        self.preq_jobs = totaljobs
        self.inq_jobs = 0
        self.completejobs = 0
        self.name = name
        
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
