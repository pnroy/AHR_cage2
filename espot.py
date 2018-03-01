from MMTK import Units
from numpy import asarray, pi, linalg, sqrt, dot

## ALL ATOM POSITIONS MUST BE IN UNITS OF : nm ##

f = Units.electrostatic_energy
charges = sqrt(f)*asarray([-0.8476, 0.4238, 0.4238, 0.4932, 0.4932, -0.9864])

def espot(atomPositions, H1pos, H2pos, Hcom):
    espot = 0.
    for i in range(3):
        H1len = atomPositions[i] - H1pos
        H2len = atomPositions[i] - H2pos
        Hcomlen = atomPositions[i] - Hcom
        espotH1 = (charges[i]*charges[3])/(sqrt(dot(H1len, H1len)))
        espotH2 = (charges[i]*charges[4])/(sqrt(dot(H2len, H2len)))
        espotHcom = (charges[i]*charges[5])/(sqrt(dot(Hcomlen, Hcomlen)))
        espot += espotH1 + espotH2 + espotHcom
    return espot
