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
ntheta = 5
nchi = 5
##############

### POTENTIAL ###
H2bl = 0.7414*Units.Ang
#################

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dz=(zmax-zmin)/nz
dtheta= pi/ntheta
dchi = 2*pi/nchi
r = H2bl/2

file_name = argv[1]
f = open(file_name, 'r')
H2Oxyzstr = list(f.readlines())
f.close()
H2Oxyz = [float(x)*Units.Ang for x in H2Oxyzstr]
H2Opositions = [H2Oxyz[x:x+3] for x in xrange(0, len(H2Oxyz), 3)] #xyz oxygen, then H1 then H2...
# xyz_out=open('cage.xyz','w')
# xyz_out.write('60'+'\n')
# xyz_out.write('a cage'+'\n')

# for a in H2Opositions:
    # xyz_out.write('O '+str(a[0])+' '+str(a[1])+' '+str(a[2])+'\n')

#####ARRAYS#####
Hcom = np.zeros(3, float)
xarr = np.fromiter([(xmin+xi*dx) for xi in xrange(nx+1)], np.float)
yarr = np.fromiter([(ymin+yi*dy) for yi in xrange(ny+1)], np.float)
zarr = np.fromiter([(zmin+zi*dz) for zi in xrange(nz+1)], np.float)
costarr = np.fromiter([(np.cos(thetai)) for thetai in xrange(ntheta+1)], np.float)
sintarr = np.fromiter([(np.sin(thetai)) for thetai in xrange(ntheta+1)], np.float)
sincarr = np.fromiter([(np.cos(chii)) for chii in xrange(nchi+1)], np.float)
coscarr = np.fromiter([(np.sin(chii)) for chii in xrange(nchi+1)], np.float)
Varr = np.zeros(((nx+1)*(ny+1)*(nz+1)), float)

Vi = 0
for xp in xrange(nx+1):
    x = xarr[xp]
    Hcom[0] = x
    for yp in xrange(ny+1):
        y = yarr[yp]
        Hcom[1] = y
        for zp in xrange(nz+1):
            z = zarr[zp]
            Hcom[2] = z
            for n in xrange(len(H2Opositions)//3):
                ljv = ljpot(H2Opositions[3*n], Hcom)
            V=0.
            for thetap in xrange(ntheta+1):
                sint = sintarr[thetap]
                cost = costarr[thetap]
                for chip in xrange(nchi):               
                    sinc = sincarr[chip]
                    cosc = coscarr[chip]
                    H2pos1 = asarray((x + r*sint*cosc), (y + r*sint*sinc), (z + r*cost))
                    H2pos2 = asarray((x - r*sint*cosc), (y - r*sint*sinc), (z - r*cost))
                    for n in xrange(len(H2Opositions)//3):
                        V+= espot(H2Opositions[3*n:3*n+3], H2pos1, H2pos2, Hcom)
            Vi += 1
            Varr[Vi-1] = (V+ljv)*dchi*dtheta/(4.*pi)

np.savetxt("optalavipot" + str(len(H2Opositions)//3), Varr, newline='\n')
