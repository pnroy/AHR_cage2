from sys import argv
from math import *
from numpy import *
from MMTK import Units
from espot import *
from ljpot import *
from scipy import special 

##################################
####### LIST OF PARAMS ###########
##################################

#### GRID ####
xmin = -1.6/10.
xmax = 1.0/10.
nx   = 5
ymin = -1.6/10.
ymax = 1.4/10.
ny   = 5
zmin = -1.7/10.
zmax = 1.7/10.
nz   = 5
ntheta = 5
nchi = 5
##############


### POTENTIAL ###
n = 3
H2bl = 0.7414*Units.Ang
#################

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dz=(zmax-zmin)/nz
dtheta = 180/ntheta
dcostheta=2./ntheta
dchi = 360/nchi
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

Hcom = zeros(3,float)
ml = [(l, m) for l in range(n+1) for m in range(-l, l+1, 1)]
mli = list(enumerate(ml))
H = zeros((n,n),complex)
# mlp = [(l, m) for l in range(n+1) for m in range(-l, l+1, 1)]


# print Hcom

potfile = open("alavipotrot" + str(len(H2Opositions)//3), "w")
for xp in range(nx+1):
    x=xmin+xp*dx
    Hcom[0] = x
#    print Hcom
    for yp in range(ny+1):
        y=ymin+yp*dy
        Hcom[1] = y
#        print Hcom
        for zp in range(nz+1):
            z=zmin+zp*dz
            Hcom[2] = z
#            print Hcom
            intgrl = 0.
            for i, ml in mli[0]:
                for j, mlp in mli[0]:
                    for theta_index in range(ntheta):
                        costheta=-1.+dcostheta*float(theta_index)
                        theta =180.*arccos(costheta)/pi
                        for chip in range(nchi+1):               
                            chi = chip*dchi
                            H2pos1 = asarray((x + r*sin(radians(theta))*cos(radians(chi))), (y + r*sin(radians(theta))*sin(radians(chi))), (z + r*cos(radians(theta))))
                            H2pos2 = asarray((x - r*sin(radians(theta))*cos(radians(chi))), (y - r*sin(radians(theta))*sin(radians(chi))), (z - r*cos(radians(theta))))
                            for n in range(len(H2Opositions)//3):
                                V = (espot(H2Opositions[3*n:3*n+3], H2pos1, H2pos2, Hcom) + ljpot(H2Opositions[3*n], Hcom))*dchi*dcostheta
                                Harm1 = special.sph_harm(mli[i][1], mli[i][0], theta, chi)
                                Harm2 = special.sph_harm(mli[j][1], mli[j][0], theta, chi)
                                print Harm1
                                C = sin(radians(theta))
                                intgrl = V*Harm1*Harm2*C*dchi*dcostheta
                                H = H[(i, j), intgrl]
            H = linalg.eig(H)            
            potfile.write(str(x)+' '+str(y)+' '+str(z)+' '+str(H/(4.*pi))+"\n")
       
potfile.close()
