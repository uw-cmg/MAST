import logging

def initialize_logger(filename="default.log"):
    logger     = logging.getLogger(filename)
    if not getattr(logger, 'handler_set', None):
        handler    = logging.FileHandler(filename)
        formatter  = logging.Formatter('%(asctime)s : %(module)15s:%(lineno)4d> : %(levelname)8s : %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.handler_set = True
    return logger
