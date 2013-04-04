from MAST.DAG.dagscheduler import DAGScheduler
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import Session
from MAST.DAG.dagscheduler import DAGParser
from MAST.DAG.dagscheduler import SessionEntry
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import JobTable
from MAST.DAG.dagscheduler import SessionTable

#from MAST.DAG.dagscheulder import Jobstatus
#from MAST.DAG.dagscheulder import Jobstatus
myds = DAGScheduler()

# Parse input file. 
# hw 2013-04-03 This is old version of input file.
depdict = DAGParser.parse('jobs_input.txt')

jt = JobTable()
st = SessionTable()
'''
    S1 = {J1, J2, J3, J4}
    Type = A  B   C   B        
    S2 =  {J1, J2}   
    Type =  A  A
'''
#myds.
#task_list = myds.schedule("../recipetest/recipe_template.txt")

