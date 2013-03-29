import os 
import errno
# root directory for experimental results
exp_root_dir = '/Users/hwkim/workspace/experiments'
# session id. Hopefully this will be assigned by DAGScheduler or something.

sid = 's001'
jid = 'j001'
spath = os.path.join(exp_root_dir,sid)
#spath = '/results'
from filemanager import FileManager
fm = FileManager()
fm.register(sid,spath)
A = 'string A'
B = 3.14159265359
fm.save_vars(sid=sid, jid=jid, gdict=globals(), varlist=['A','B'])
del A, B
#load variables to workspace
fm.load_variables_to_ws(sid, jid, globals())

#load variable dictionary
vardict = fm.load_variables(sid, jid)
