from sys import argv
import numpy as np
from MMTK import Units
from espot import *
from ljpot import *

##################################
####### LIST OF PARAMS ###########
##################################

#### GRID ####
xmin = -5.6*Units.Ang
xmax = 5.0*Units.Ang
nx   = 5
ymin = -5.6*Units.Ang
ymax = 5.4*Units.Ang
ny   = 5
zmin = -5.7*Units.Ang
zmax = 5.7*Units.Ang
nz   = 5
ncostheta = 5
nchi = 5
##############

### POTENTIAL ###
H2bl = 0.7414*Units.Ang
#################

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dz=(zmax-zmin)/nz
dcostheta=2./ncostheta
dchi = 2*pi/nchi
r = H2bl/2.

file_name = argv[1]
f = open(file_name, 'r')
H2Oxyzstr = list(f.readlines())
f.close()
H2Oxyz = [float(x)*0.1 for x in H2Oxyzstr]
H2Opositions = [H2Oxyz[x:x+3] for x in xrange(0, len(H2Oxyz), 3)] #xyz oxygen, then H1 then H2...
# xyz_out=open('cage.xyz','w')
# xyz_out.write('60'+'\n')
# xyz_out.write('a cage'+'\n')

# for a in H2Opositions:
    # xyz_out.write('O '+str(a[0])+' '+str(a[1])+' '+str(a[2])+'\n')

#####ARRAYS#####
Hcom = np.zeros(3,float)




potfile = open("alavipot" + str(len(H2Opositions)//3), "w")
for xp in range(nx+1):
    x=xmin+xp*dx
    Hcom[0] = x
    for yp in range(ny+1):
        y=ymin+yp*dy
        Hcom[1] = y
        for zp in range(nz+1):
            z=zmin+zp*dz
            Hcom[2] = z
            V=0.
            for costheta_index in range(ncostheta+1):
                cost=-1.+dcostheta*float(costheta_index)
                theta = np.arccos(cost)
                sint = np.sin(theta)
                for chip in range(nchi):               
                    chi = chip*dchi
                    sinc = np.sin(chi)
                    cosc = np.cos(chi)
                    H2pos1 = asarray((x + r*sint*cosc), (y + r*sint*sinc), (z + r*cost))
                    H2pos2 = asarray((x - r*sint*cosc), (y - r*sint*sinc), (z - r*cost))
                    for n in range(len(H2Opositions)//3):
                        V+= espot(H2Opositions[3*n:3*n+3], H2pos1, H2pos2, Hcom) + ljpot(H2Opositions[3*n], Hcom)
            potfile.write(str(x)+' '+str(y)+' '+str(z)+' '+str(V*dchi*dcostheta/(4.*pi))+"\n")

potfile.close()
