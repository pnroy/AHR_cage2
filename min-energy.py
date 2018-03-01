from sys import argv
from numpy import argmin

f = argv[1]
ifile = open(f, 'r')
minval=0.
for i in range(216):
    energies = ifile.readline()
    print energies
    if (energies[0]==''):
        exit
    energies=float(energies.split()[3])
    if (energies < minval):
        minval = energies
        minind=i


print minind, minval
ifile.close()
