from sys import argv
from scipy import integrate, special
import numpy as np
from MMTK import Units
from numpy import asarray, pi, linalg, sqrt, dot

## Certain sections may be unnecessary currently, still tweaking ##

##################################
####### LIST OF PARAMS ###########
##################################

### POTENTIAL ###
n = 3
H2bl = 0.7414*Units.Ang
r = H2bl/2.
#################

#### POTENTIAL FUNCTION DEFINITIONS ####
## ALL ATOM POSITIONS MUST BE IN UNITS OF : nm ##

## Electrostatic contribution - 9 interactions between each H2 and H20.
## Takes in position 
def espot(atomPositions, H1pos, H2pos, Hcom):
    f = Units.electrostatic_energy
    charges = sqrt(f)*asarray([-0.8476, 0.4238, 0.4238, 0.4932, 0.4932, -0.9864])
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


def ljpot(atomPositions, Hcom):
    sig = 3.101339711*Units.Ang
    eps = 0.4306
    HcomM = atomPositions - Hcom 
    r2=(dot(HcomM,HcomM))
    r6=r2*r2*r2
    r12=r6*r6
    V = 4.*eps*((sig)**12./r12 - (sig)**6./r6)
    return V

## Reads in positions from file (required argument) ##
file_name = argv[1]
f = open(file_name, 'r')
H2Oxyzstr = list(f.readlines())
f.close()
H2Oxyz = [float(x)*0.1 for x in H2Oxyzstr]
H2Opositions = [H2Oxyz[x:x+3] for x in xrange(0, len(H2Oxyz), 3)] #xyz oxygen, then H1 then H2...
# xyz_out=open('cage.xyz','w')
# xyz_out.write('60'+'\n')
# xyz_out.write('a cage'+'\n')
## setup basis size
sizej=2;
jmax=2*(sizej-1);

## grid size

size_theta=2*jmax+5
size_phi=2*(2*jmax+7)

theta_points, theta_weights =np.polynomial.legendre.leggauss(size_theta)

outfile=open('grid_theta','w')

index_point=0
for point in theta_points:
	outfile.write(str(point)+' '+str(theta_weights[index_point])+'\n')
	index_point+=1

phi_points=np.zeros(size_phi,'float')

delta_phi=2.*np.pi/float(size_phi)

for i in range(size_phi):
	phi_points[i]=float(i)*delta_phi

index_point=0
for point in phi_points:
	outfile.write(str(point)+'\n')
        index_point+=1

for j in range(0,jmax+1,2):
	for m in range(-j,j+1,1):
		for jp in range(0,jmax+1,2):
			for mp in range(-jp,jp+1,1):
				test_ortho=complex(0.,0.)
				for k in range(len(theta_points)):
					theta=np.arccos(theta_points[k])
					for phi in phi_points:
						test_ortho+=theta_weights[k]*delta_phi*np.conj(special.sph_harm(m, j, phi, theta))*special.sph_harm(mp, jp, phi,theta)
				print j,m,jp,mp,test_ortho
