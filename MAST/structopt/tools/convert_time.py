import math

def convert_time(t):
    """Function to convert time in seconds to readable string"""
    tdays = math.floor(t/(3600*24.0))
    tleft = t-tdays*3600*24
    thours = math.floor(tleft/3600)
    tleft -=thours*3600
    tmins = math.floor(tleft/60)
    tleft -=tmins*60
    tseconds = tleft
    str = "%02d:%02d:%02d:%02.4f" % (tdays, thours, tmins, tseconds)
    return str

