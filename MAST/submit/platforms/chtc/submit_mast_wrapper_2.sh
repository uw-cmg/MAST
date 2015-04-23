#!/bin/sh
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o mastdeps.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/mastdeps.tar.gz
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o install_mast_on_epd.sh http://proxy.chtc.wisc.edu/SQUID/mayeshiba/install_mast_on_epd.sh
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o MAST_workdirs_2.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/MAST_workdirs_2.tar.gz
curl -H "Pragma:" --fail --retry 30 --retry-delay 6 -o vsquid_pps.tar.gz http://proxy.chtc.wisc.edu/SQUID/mayeshiba/vsquid_pps.tar.gz

chmod u+x install_mast_on_epd.sh
./install_mast_on_epd.sh
tar -xzvf MAST_workdirs_2.tar.gz
tar -xzvf vsquid_pps.tar.gz

export PATH=$TEMP/epd/bin:$PATH
export MAST_SCRATCH=$TEMP/MAST_2/SCRATCH
export MAST_ARCHIVE=$TEMP/MAST_2/ARCHIVE
export MAST_CONTROL=$TEMP/MAST_2/CONTROL
export MAST_PLATFORM=chtc_write_only
export VASP_PSP_DIR=$TEMP/vasp_pps 


rm $TEMP/MAST_2/SCRATCH/mast.write_files.lock
./comparescratch_2.py -m save
mast -m monitoronly > mastoutput
sleep 3600
./comparescratch_2.py -m compare
./comparescratch_2.py -m zip
rm -r $MAST_CONTROL/statusfiles
tar -czvf back_2.tar.gz -T tarthis
rm vsquid_pps.tar.gz
rm mastdeps.tar.gz
rm MAST_workdirs_2.tar.gz
rm -r vasp_pps
rm -r epd
rm -r MAST
exit
