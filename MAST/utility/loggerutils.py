import logging

def initialize_logger(filename="default.log", formatstr=""):
    logger     = logging.getLogger(filename)
    if not getattr(logger, 'handler_set', None):
        handler    = logging.FileHandler(filename)
        if formatstr == "":
            #formatter  = logging.Formatter('%(asctime)s : %(module)15s:%(lineno)4d> : %(levelname)8s : %(message)s')
            formatter  = logging.Formatter('%(asctime)s : %(module)10s: %(levelname)8s : %(message)s')
        else:
            formatter = logging.Formatter(formatstr)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.handler_set = True
    return logger

def initialize_short_logger(filename="default.log"):
    format = "%(module)10s: %(message)s"
    #format = "%(asctime)s %(levelname)8s : %(message)s"
    return initialize_logger(filename, format)
