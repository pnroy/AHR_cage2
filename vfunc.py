from sys import argv
import numpy as np
from MMTK import Units

### GRID ###
xmin = -5.6*Units.Ang
xmax = 5.0*Units.Ang
nx   = 1000
ymin = -5.6*Units.Ang
ymax = 5.4*Units.Ang
ny   = 5
zmin = -5.7*Units.Ang
zmax = 5.7*Units.Ang
nz   = 5
ntheta = 5
nchi = 5
dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dz=(zmax-zmin)/nz
dtheta= np.pi/ntheta
dchi = 2*(np.pi)/nchi

### Constants ###
H2bl = 0.7414*Units.Ang
sig = 3.101339711*Units.Ang
eps = 0.4306
f = Units.electrostatic_energy
qqlst = np.sqrt(f)*[-0.4180, 0.2090, 0.2090, 0.8361, -0.4180, -0.4180]
r = H2bl/2

### Computing H2O positions ###
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

### LJ potential function ###
def ljpot(atomPos, Hcom):
        rVec = atomPos - Hcom
        r2 = np.dot(rVec, rVec)
        r6 = r2*r2*r2
        r12 = r6*r6
        ljpot = 4.*eps*((sig)**12./r12 - (sig)**6./r6)
        return ljpot

### Electrostatic interactions function ###
def espot(atomPos, H1pos, H2pos, Hcom):
        espot = 0.
        for i in range(3):
            H1R = atomPos[i] - H1pos
            H2R = atomPos[i] - H2pos
            HcomR = atomPos[i] - Hcom
            H1Len = np.sqrt(np.dot(H1R,H1R))
            H2Len = np.sqrt(np.dot(H2R,H2R))
            HcomLen = np.sqrt(np.dot(HcomR,HcomR))
            espotH1 = qqlst[i]/H1Len
            espotH2 = qqlst[i]/H2Len
            espotHcom = qqlst[i+3]/HcomLen
            espot += espotH1 + espotH2 + espotHcom
        return espot

### Overall potential function ###
def vfunc(t, p, x, y, z):
    V = 0.
    ### LJ ###
    Hcom = np.asarray([x, y, z])
    for n in xrange(len(H2Opositions)//3):
        V += ljpot(H2Opositions[3*n], Hcom)
    ### ES ###
    sint = np.sin(t)
    cost = np.cos(t)
    sinp = np.sin(p)
    cosp = np.cos(p)
    H2pos1 = np.asarray((x + r*sint*cosp), (y + r*sint*sinp), (z + r*cost))
    H2pos2 = np.asarray((x - r*sint*cosp), (y - r*sint*sinp), (z - r*cost))  
    for n in xrange(len(H2Opositions)//3):
        V += espot(H2Opositions[3*n:3*n+3], H2pos1, H2pos2, Hcom)
    
    return V

### Quick testing ###
xarr = np.fromiter([(xmin+xi*dx) for xi in xrange(nx+1)], np.float)
Varr = np.zeros(nx+1)
Vi = 0
for xi in xrange(nx + 1):
    V = 0.
    x = xarr[xi]
    V = vfunc(0., 0., x, 0., 0.)
    Vi += 1
    Varr[Vi-1] = V
    
np.savetxt("vfunctest" + str(len(H2Opositions)//3), Varr, newline='\n')

