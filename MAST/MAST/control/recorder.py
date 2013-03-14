############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################

#20130128 TTM created to house common control functions; fcns are identical
#             and should probably be reevaluated for organization

import os
import time
#import smtplib
from MAST.mastfile.mastfile import MASTFile

class Recorder:
    """Recorder class allows recording of output.
    
    ???There must be a better way to do this.

    Attributes:
        controlpath
    """
    def __init__(self):
        self.controlpath = "" #<str>: path to control folder
        self.controlpath = os.path.join(os.path.expanduser(os.getenv("SCRIPTPATH")),"controlfolder")
        #import MAST.utility.globals

    def record_output(self, mypath, outputmsg):
        """Records an output message to output.txt.
            
        Args:
            mypath <str>: path in which to create or write output.txt
            outputmsg <str>: message to write
        """
        fullpath=""
        #fullpath = os.path.join(mypath,time.strftime('%Y%m%d%H%M%S')+'_output.txt')
        fullpath = os.path.join(mypath,'output.txt')
        ofile = MASTFile(fullpath)
        ofile.data.append(time.asctime())
        ofile.data.append(' ' + outputmsg + '\n') #TTM 20121024 add space
        ofile.to_file(fullpath)
        return
    def send_message(self, emsubject, emmessage):
        """Sends an email message.
        
        Args:
            emsubject <str>: subject line of email
            emmessage <str>: message of email.
        """
        #TTM 052212 created
        sess=""
        msg=""
        msg = email.message_from_string(emmessage)
        msg['Subject'] = emsubject
        msg['From'] = 'MAterials Simulation Toolbox (Automatam)'
        msg['To'] = email_path #global variable
        sess=smtplib.SMTP("localhost")
        sess.sendmail(email_path, email_path, msg.as_string())
        sess.quit()
        return True        
