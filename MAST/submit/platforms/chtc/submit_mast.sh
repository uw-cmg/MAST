################################
universe = vanilla
executable = //home/tmayeshi/bin/submit_mast_wrapper.sh
transfer_input_files =  //home/tmayeshi/bin/comparescratch.py, //home/tmayeshi/canopy_local.tar.gz, //home/tmayeshi/MAST_workdirs.tar.gz, //home/tmayeshi/vasp_pps.tar.gz
should_transfer_files = yes
when_to_transfer_output = on_exit
request_cpus = 1
request_memory = 1000
request_disk = 1000000

#Requirements = CAN_RUN_WHOLE_MACHINE

log = log
output = output
error  = error
getenv = true
#+RequiresWholeMachine=true
#+AccountingGroup = MSE_Morgan
+WantFlocking = True
+WantGlideIn = True
queue
################################
