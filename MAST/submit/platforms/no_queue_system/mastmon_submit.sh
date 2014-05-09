#!/bin/bash
nice -n 19 python runmast.py >> $MAST_CONTROL/mastoutput 2> $MAST_CONTROL/errormast
