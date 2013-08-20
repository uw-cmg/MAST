#!/usr/bin/env python
from MAST.DAG.mastmon import MASTmon

mastmon = MASTmon()
mastmon.run(1,1)
#if int(mastopt.verbose) == 1:
#    mastmon.run(1,1) #TTM run for only 1 dagscheduler iteration. Use crontab to check periodically.
#else:
#    mastmon.run(1,0)

