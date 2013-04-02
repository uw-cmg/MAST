'''All files for a session will be stored in a directory and subdirectories of it.
    Ex1)
    experiment1/
        ./job1/
        ./job2/
        ./job3/

    Ex2) <- the simplest way. I implemented this

    experiment2/
        ./job1.pkl
        ./job2.pkl
        ./job3.pkl
    Ex3)
    experiment3/
        ./jid_01.pkl  # I think some jobs might be done  many time in the loop.
        ./jid_02.pkl  # jid maybe different even though job name is same. 
        ./jid_03.pkl  # if number of iteration is big then the intermediate files blow up
'''

import os 
import errno
from picklemanager import PickleManager

class FileManager(PickleManager):
    '''FILEMANAGER will load correct pickle files by session id (sid) and job id (jid)'''

    def __init__(self):
        self.dirdict = {}
        PickleManager.__init__(self)

    def register(self, sid, session_dir):
        if self.makedir(session_dir):
            self.dirdict[sid] = session_dir
        else:
            raise Exception("Failed to make a session directory. \nDetail: %s" % ex.strerror)

    def makedir(self, path):
        '''MAKEDIR makes a directory (path) or check that the 
            directory exist, then return True.
            If it fails to make a directory, then returns False.
        '''
        try:
            os.makedirs(path)
            print spath
        except OSError as ex:
            if ex.errno == errno.EEXIST:
                print 'Results will be overwritten in %s' % path
            else:
                print ex.strerror
                return False
        return True

    def save_vars(self, sid, jid, gdict, varlist=None):
        '''SAVE_VAR saves variables in a pickle file in a session directory'''
        fURL = self.getAbsURL(sid,jid)
        self.save_variables(gdict=gdict, varlist=varlist, filename=fURL)

    def load_variables(self, sid, jid):
        '''LOAD_VARIABLES returns a dictionary of variables.
            sid : session id
            jid : job id
        '''
        fURL = self.getAbsURL(sid,jid)
        return PickleManager.load_variables(self, filename=fURL)
    def load_variables_to_ws(self, sid, jid, gdict):
        fURL = self.getAbsURL(sid,jid)
        return PickleManager.load_variables_to_ws(self, gdict=gdict, filename=fURL)
        
    def getAbsURL(self, sid, jid):
        '''GETFILEPATH returns full url of the picke file of (sid,jid)
        '''
        fURL = None
        if sid in self.dirdict:
            fURL = os.path.join(self.dirdict[sid], jid+'.pkl')   
        return fURL











