################################
universe = vanilla
executable = //home/tmayeshi/bin/wrapper_exec
arguments = ?mast_exec? ?mast_ppn?
transfer_input_files = ./,?mast_exec?,//home/tmayeshi/openmpi
should_transfer_files = yes
when_to_transfer_output = on_exit
request_cpus = ?mast_ppn?
request_memory = 1000
request_disk = 1000

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
