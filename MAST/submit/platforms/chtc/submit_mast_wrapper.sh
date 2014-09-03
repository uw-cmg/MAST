#!/bin/sh
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o canopy_local.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/canopy_local.tar.gz
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o MAST_workdirs.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/MAST_workdirs.tar.gz
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o vsquid_pps.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/vsquid_pps.tar.gz

tar -xzvf canopy_local.tar.gz
tar -xzvf MAST_workdirs.tar.gz
tar -xzvf vsquid_pps.tar.gz

export PATH=$TEMP/canopy_local/appdata/canopy-1.4.1.1975.rh5-x86_64/bin:$PATH
export MAST_SCRATCH=$TEMP/MAST/SCRATCH
export MAST_ARCHIVE=$TEMP/MAST/ARCHIVE
export MAST_CONTROL=$TEMP/MAST/CONTROL
export MAST_PLATFORM=no_queue_system
export MAST_RECIPE_PATH=$TEMP/MAST/recipe_templates
export VASP_PSP_DIR=$TEMP/vasp_pps 

./comparescratch.py -m save
mast -m monitoronly
sleep 100
./comparescratch.py -m compare
./comparescratch.py -m zip
rm -r $MAST_CONTROL/statusfiles
tar -czvf back.tar.gz -T tarthis
rm -r vasp_pps
rm -r canopy_local
rm -r MAST
exit
