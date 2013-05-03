from jobentry import JobEntry
from dagutil import *

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
            
