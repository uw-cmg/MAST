from MAST.DAG.mastmon import MASTmon
import os
import glob

mastmon = MASTmon()
mastmon.run()

#            ex) mastmon.run()  # run MASTmon forever as a real daemon. By default interval is 10 sec
#            ex) mastmon.run(interval=30) # run MASTmon forever as a real daemon. By default interval is 30 sec
#            ex) mastmon.run(niter=1) # run MASTmon one iteration for crontab user. By default interval is 10 sec
#            ex) mastmon.run(niter=20,stopcond='NOSESSION') # run MASTmon for 20 iterations. And stop it all sessions are done.
