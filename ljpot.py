#Lennard-Jones Potential
from numpy import sqrt,dot
from MMTK import Units

## ATOM POSITIONS MUST BE IN : nm ##

sig = 3.101339711*Units.Ang
eps = 0.4306

def ljpot(atomPositions, Hcom):
    HcomM = atomPositions - Hcom 
    r2=(dot(HcomM,HcomM))
    r6=r2*r2*r2
    r12=r6*r6
    V = 4.*eps*((sig)**12./r12 - (sig)**6./r6)
    return V
