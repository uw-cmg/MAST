#!/bin/bash
tar -xzvf mastdeps.tar.gz
cd mastdeps
./epd_free-7.3-2-rh5-x86_64.sh -b -p $TMP/epd
export PATH=$TMP/epd/bin:$PATH
while read line; do echo $line; sname=`echo $line | awk -F. '{print $1}'`; echo $sname; tar -xzvf $line; cd $sname*; python setup.py install; cd $TMP/mastdeps; done < deplist
