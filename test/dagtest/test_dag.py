from MAST.DAG.dagscheduler import DAGScheduler
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import DAGParser
from MAST.DAG.dagscheduler import SessionEntry
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import JobTable
from MAST.DAG.dagscheduler import SessionTable
from MAST.DAG.dagscheduler import set2str
from sys import exit

A = set()
A.add('dog')
A.add('cat')
A.add('pig')
A.add('rabbit')
out = set2str(A)
print out


#from MAST.DAG.dagscheulder import Jobstatus
#from MAST.DAG.dagscheulder import Jobstatus
myds = DAGScheduler()

# Parse input file. 
# hw 2013-04-03 This is old version of input file.
depdict = DAGParser.parse('jobs_input.txt')

jt = JobTable()
st = SessionTable()

j1 = JobEntry(1,1,'job1','./','./','A')
j2 = JobEntry(1,2,'job2','./','./','B')
j3 = JobEntry(1,3,'job3','./','./','C')
j4 = JobEntry(1,4,'job4','./','./','D')

jt.addjob(j1)
jt.addjob(j2)
jt.addjob(j3)
jt.addjob(j4)

# jid confliction raises exception
try:
    jt.addjob(j1)
except Exception as ex:
    print ex.message

#try to delete a job which is not in jobtable
try:
    jt.deljob(5)
except Exception as ex:
    print ex.message

print 'Job table print test'
print jt

# clear jobs by jid    
for id in range(4):
    jt.deljob(id+1)

print 'Session table print test'
s1 = SessionEntry(1,4)
s2 = SessionEntry(2,3)
st.addsession(s1)
st.addsession(s2)

print 'After deletion'
print jt

# To print out table
print st


    
'''
    S1 = {J1, J2, J3, J4}
    Type = A  B   C   B        
    S2 =  {J1, J2}   
    Type =  A  A
'''
#myds.
#task_list = myds.schedule("../recipetest/recipe_template.txt")

#TODO LIST for DAG Scheduler
# Make a class manage tables and submit jobs via PBS for something
# Make a class or module to communicate with MASTMon or Daemon.
