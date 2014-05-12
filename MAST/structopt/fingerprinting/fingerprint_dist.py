import numpy

def fingerprint_dist(fingerprint1, fingerprint2):
    """Function to calculate the cosine distance between two fingerprint functions
    """
    f1=numpy.array(fingerprint1)
    f1mag=numpy.linalg.norm(f1)
    f2=numpy.array(fingerprint2)
    f2mag=numpy.linalg.norm(f2)
    dist=0.5*(1-numpy.dot(f1,f2)/(f1mag*f2mag))
    return dist

