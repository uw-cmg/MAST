################################
universe = vanilla
executable = submit_mast_wrapper_2.sh
transfer_input_files = comparescratch_2.py
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
#+WantFlocking = True
#+WantGlideIn = True
queue
################################
