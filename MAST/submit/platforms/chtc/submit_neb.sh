################################
universe = vanilla
executable = //home/mayeshiba/bin/wrapper_exec_vsquid_checkpoint
arguments = 16
transfer_input_files = ./
transfer_output_files = vasprun.xml,00,01,02
should_transfer_files = yes
when_to_transfer_output = ON_EXIT_OR_EVICT
+is_resumable = true
request_cpus = 16
request_memory = 20000
request_disk = 10000000

log = log
output = output
error  = error
getenv = true
#+AccountingGroup = MSE_Morgan
#+WantFlocking = True
#+WantGlideIn = True
queue
################################
