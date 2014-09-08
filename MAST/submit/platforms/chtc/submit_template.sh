################################
universe = vanilla
executable = //home/mayeshiba/bin/wrapper_exec_?mast_exec?
arguments = ?mast_ppn?
transfer_input_files = ./
should_transfer_files = yes
when_to_transfer_output = on_exit
request_cpus = ?mast_ppn?
request_memory = 5000
request_disk = 1000000

Requirements = CAN_RUN_WHOLE_MACHINE

log = log
output = output
error  = error
getenv = true
+RequiresWholeMachine=true
#+AccountingGroup = MSE_Morgan
+WantFlocking = True
+WantGlideIn = True
queue
################################
