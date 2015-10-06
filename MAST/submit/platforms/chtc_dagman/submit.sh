################################
universe = vanilla
executable = ./wrapper_exec_vsquid_checkpoint
arguments = 16
transfer_input_files = ./
should_transfer_files = yes
when_to_transfer_output = ON_EXIT_OR_EVICT
+is_resumable = true
request_cpus = 16
request_memory = 100000
request_disk = 10000000

#Requirements = CAN_RUN_WHOLE_MACHINE

log = log
output = output
error  = error
getenv = true
#+RequiresWholeMachine=true
#+AccountingGroup = MSE_Morgan
#+WantFlocking = True
#+WantGlideIn = True
queue
################################
