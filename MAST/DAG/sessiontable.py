from dagutil import *
from sessionentry import SessionEntry

class SessionTable(object):
    def __init__(self):
        self.sessions = {}
        self._nextsid = 1
        self._sname2sid ={}
        
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
        self._sname2sid[session.sname] = session.sid
        
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
