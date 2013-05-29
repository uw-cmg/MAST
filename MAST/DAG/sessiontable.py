from dagutil import *
from sessionentry import SessionEntry

class SessionTable(object):
    """SessionTable has multiple SessionEntry objects. Those are for keeping information of each session
        Each session is independant each other.
    """
    def __init__(self):
        self.sessions = {}
        self._nextsid = 1
        self._sname2sid ={}
        
    def completejobs(self, sid, njobs):
        """ Update SessionEntry sid to set njobs are complete (Complete). """
        self.sessions[sid].set_jobs_to_complete(njobs)
        
    def submitjobs(self, sid, njobs):
        """ Update SessionEntry sid to set njobs are submitted (PreQ). """
        self.sessions[sid].set_jobs_to_submitted(njobs)

    def runjobs(self, sid, njobs):
        """ Update SessionEntry sid to set njobs are running in queue (InQ). """
        self.submitjobs(sid, njobs)
        
    def addsession(self, session):
        """ Add sessions to SessionTable."""
        if session.sid in self.sessions:
            raise Exception("SESSION ID (sid=%d) CONFLICTED" % session.sid)
        self.sessions[session.sid] = session
        self._sname2sid[session.name] = session.sid
        
    def delsession(self, sid):
        """ Delete Session sid from SessionTable."""
        if sid not in self.sessions:
            raise Exception("SESSION ID (sid=%d) DOESN'T EXIST" % sid)
        session = self.sessions.pop(sid)
        del self._sname2sid[session.name]
        
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
    
    def get_sid(self, sname):
        return self._sname2sid[sname]
    
    def get_sname(self,sid):
        return self.sessions[sid].name
    
    def get_new_sid(self):
        """ get_sid returns unique session id in a session table object. """
        while (self._nextsid in self.sessions):
            self._nextsid = (self._nextsid ) % MAXSID + 1
        return self._nextsid

    def add_newjob(self, sid, njobs, status = JOB.PreQ):
        """ Add jobs to SessionEntry sid in SessionTable. """
        self.sessions[sid].add_newjob(njobs,status)
