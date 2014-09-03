################################
universe = vanilla
executable = //home/mayeshiba/bin/submit_mast_wrapper.sh
transfer_input_files =  //home/mayeshiba/bin/comparescratch.py
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
