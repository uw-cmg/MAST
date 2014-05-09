def calc_dist(atom1,atom2):
    """Function to calculates the x, y, z, and total distance between atoms
    Inputs:
    	atom1 = ASE class atom object
    	atom2 = ASE class atom object
    Outputs:
    	[d,x,y,z] = [total distance, x-direction distance, y-direction distance, z-direction distance]
    	*Note: all distances are in Angstroms
    """
    pos1=atom1.position
    pos2=atom2.position
    x=abs(pos1[0]-pos2[0])
    y=abs(pos1[1]-pos2[1])
    z=abs(pos1[2]-pos2[2])
    d=(x**2+y**2+z**2)**0.5
    return [d,x,y,z]
