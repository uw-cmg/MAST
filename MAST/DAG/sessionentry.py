from dagutil import *

class SessionEntry(object):

    def __init__(self, sid, totaljobs, name=None):
        self.sid =  sid
        self.totaljobs = totaljobs
        self.preq_jobs = totaljobs
        self.inq_jobs = 0
        self.completejobs = 0
        self.name = name
        self.preq_list=list()
        self.inq_list=list()
        self.complete_list=list()
        self.total_list=list()
    
    def count_lists(self):
        self.inq_jobs = len(self.inq_list)
        self.preq_jobs = len(self.preq_list)
        self.completejobs = len(self.complete_list)
        self.totaljobs = len(self.total_list)
        self._print_lists()

    def _print_lists(self):
        #pass
        print "INQ:",self.inq_list
        print "PREQ:",self.preq_list
        print "COMP:",self.complete_list
        print "TOT:",self.total_list

    def set_jobs_to_submitted(self, njobs):
        if njobs in self.preq_list:
            self.preq_list.remove(njobs)
        if not njobs in self.inq_list:
            self.inq_list.append(njobs)
        self.count_lists()

    def set_jobs_to_complete(self, njobs):
        if njobs in self.inq_list:
            self.inq_list.remove(njobs)
        if not njobs in self.complete_list:
            self.complete_list.append(njobs)
        self.count_lists()

    def add_newjob(self, njobs, status = JOB.PreQ):
        if njobs in self.total_list:
            return
        self.total_list.append(njobs)
        if status is JOB.PreQ:
            self.preq_list.append(njobs)
        elif status is JOB.InQ:
            self.inq_list.append(njobs)
        elif status is JOB.Complete:
            self.complete_list.append(njobs)
        self.count_lists()
        
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
