from scipy import integrate, special 
import numpy as np
from numpy.polynomial import chebyshev, legendre
from numpy import ndarray
from math import pi

def func(x):
    return np.exp(-5*x);

def func2(x):
    return 1/(np.sqrt(1-x**2));

def func3(x):
    return x;

val = integrate.quadrature(func, -1, 1)
print val

xarr, yarr = chebyshev.chebgauss(5)
intval = 0

nsteps = xarr.shape[0]

for n in np.arange(nsteps):
    intval += yarr[n]*func3(xarr[n])
print intval

degree = 1
theta_nodes, theta_weights = legendre.leggauss((degree+1))
phi_nodes, phi_weights = chebyshev.chebgauss(2*degree-1)


for n in range(len(xarr)):
        
    special.sph_harm(1, 1, pi*phi + pi, cos(theta))
