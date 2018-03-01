#!/usr/bin/env python
# -*- coding: utf-8 -*-

# q-TIP4P water potential
# Kevin Bishop

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'
from libc.math cimport sin,sqrt,pow,acos,rint
cimport numpy as N
import numpy as N

cdef class HarmonicAngleTerm(EnergyTerm):

    cdef int a1,a2,a3, periodic
    cdef double theta_eq,k_theta,cellSize[3]
    cdef scale

    def __init__(self, universe, atom_indices, theta_eq, k_theta, beg_end):
        cdef ArrayType ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_angle", ("harmonic_angle",))
        self.eval_func = <void *>HarmonicAngleTerm.evaluate
        #Middle atom is a2!
        self.a1 = atom_indices[0]
        self.a2 = atom_indices[1]
        self.a3 = atom_indices[2]
        self.theta_eq = theta_eq
        self.k_theta = k_theta
        self.periodic = universe.is_periodic
        self.scale = 1.
        if beg_end:
            self.scale = 0.5
        if self.periodic:
            self.cellSize[0] = universe.cellParameters()[0]
            self.cellSize[1] = universe.cellParameters()[1]
            self.cellSize[2] = universe.cellParameters()[2]

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef int i_loop
        cdef double a1_pos[3],a2_pos[3],a3_pos[3],r12[3],r32[3],d1[3],d2[3],d3[3]
        cdef double lr12 = 0.0, lr32 = 0.0
        cdef double r12_dot_r32 = 0.0
        cdef double cos_theta, sin_theta, theta, deriv
        #a2 is middle atom
        coordinates = <vector3 *>input.coordinates.data

        for i_loop in range(3):
            a1_pos[i_loop] = coordinates[self.a1][i_loop]
            a2_pos[i_loop] = coordinates[self.a2][i_loop]
            a3_pos[i_loop] = coordinates[self.a3][i_loop]
            r12[i_loop] = a1_pos[i_loop] - a2_pos[i_loop]
            r32[i_loop] = a3_pos[i_loop] - a2_pos[i_loop]
            if self.periodic:
                r12[i_loop] -= rint(r12[i_loop]/self.cellSize[i_loop])*self.cellSize[i_loop]
                r32[i_loop] -= rint(r32[i_loop]/self.cellSize[i_loop])*self.cellSize[i_loop]
            r12_dot_r32 += r12[i_loop]*r32[i_loop]
            lr12 += pow(r12[i_loop],2)
            lr32 += pow(r32[i_loop],2)
        lr12 = sqrt(lr12)
        lr32 = sqrt(lr32)

        cos_theta = r12_dot_r32 / (lr12*lr32)
        if cos_theta > 1.: cos_theta = 1.
        if cos_theta < -1.: cos_theta=-1.
        theta = acos(cos_theta)
        dtheta = theta - self.theta_eq

        #Add energy to MMTK
        #Scaling is performed on end beads for PIGS
        energy.energy_terms[self.index] = self.scale*0.5*self.k_theta*(pow(dtheta,2.))

        #Add gradients
        if energy.gradients != NULL:
            #Determine dV wrt theta
            sin_theta = sqrt(1. - pow(cos_theta,2.))
            deriv = -1.*self.k_theta*dtheta/sin_theta * self.scale   #PIGS scaling
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

            #Determine dtheta wrt cartesian component
            for i_loop in range(3):
                d1[i_loop] = (r32[i_loop]/(lr12*lr32)) - (cos_theta*r12[i_loop]/(pow(lr12,2.0)))
                d3[i_loop] = (r12[i_loop]/(lr12*lr32)) - (cos_theta*r32[i_loop]/(pow(lr32,2.0)))
                d2[i_loop] = -d1[i_loop]-d3[i_loop]
            for i_loop in range(3):
                gradients[self.a1][i_loop] += d1[i_loop]*deriv
                gradients[self.a2][i_loop] += d2[i_loop]*deriv
                gradients[self.a3][i_loop] += d3[i_loop]*deriv
  
cdef class HarmonicBondTerm(EnergyTerm):

    cdef int a1,a2,a3,periodic
    cdef double r_eq, k_r, cellSize[3]
    cdef scale

    def __init__(self, universe, atom_indices, r_eq, k_r, beg_end):
        cdef ArrayType ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_bond", ("harmonic_bond",))
        self.eval_func = <void *>HarmonicBondTerm.evaluate
        self.a1 = atom_indices[0]
        self.a2 = atom_indices[1]
        self.r_eq = r_eq
        self.k_r = k_r
        self.periodic = universe.is_periodic
        self.scale = 1.
        if beg_end:
            self.scale = 0.5
        if self.periodic:
            self.cellSize[0] = universe.cellParameters()[0]
            self.cellSize[1] = universe.cellParameters()[1]
            self.cellSize[2] = universe.cellParameters()[2]

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef int i_loop
        cdef double a1_pos[3],a2_pos[3],r12[3],d1[3],d2[3]
        cdef double lr12 = 0.0, lr32 = 0.0, r12_dot_r32 = 0.0, deriv
        coordinates = <vector3 *>input.coordinates.data

        for i_loop in range(3):
            a1_pos[i_loop] = coordinates[self.a1][i_loop]
            a2_pos[i_loop] = coordinates[self.a2][i_loop]
            r12[i_loop] = a2_pos[i_loop] - a1_pos[i_loop]
            if self.periodic:
                r12[i_loop] -= rint(r12[i_loop]/self.cellSize[i_loop])*self.cellSize[i_loop]
            lr12 += pow(r12[i_loop],2)
        lr12 = sqrt(lr12)
        dr = lr12 - self.r_eq

        #Add energy to MMTK, PIGS
        energy.energy_terms[self.index] = self.scale*0.5*self.k_r*(pow(dr,2.))

        #Add gradients
        if energy.gradients != NULL:
            deriv = -1.*self.k_r*dr * self.scale #PIGS
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

            for i_loop in range(3):
                d1[i_loop] = r12[i_loop] / lr12
                d2[i_loop] = - d1[i_loop]

            for i_loop in range(3):
                gradients[self.a1][i_loop] += d1[i_loop]*deriv
                gradients[self.a2][i_loop] += d2[i_loop]*deriv

cdef class QuarticBondTerm(EnergyTerm):

    cdef int a1,a2,periodic
    cdef double D_r,c1,c2,c3,r_eq,cellSize[3]
    cdef scale

    def __init__(self, universe, atom_indices, D_r, c1, c2, c3, r_eq, beg_end):
        cdef ArrayType ref_array
        EnergyTerm.__init__(self, universe,
                            "quartic_bond", ("quartic_bond",))
        self.eval_func = <void *>QuarticBondTerm.evaluate
        self.a1 = atom_indices[0]
        self.a2 = atom_indices[1]
        self.D_r = D_r
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.r_eq = r_eq
        self.periodic = universe.is_periodic
        self.scale = 1.
        if beg_end:
            self.scale = 0.5
        if self.periodic:
            self.cellSize[0] = universe.cellParameters()[0]
            self.cellSize[1] = universe.cellParameters()[1]
            self.cellSize[2] = universe.cellParameters()[2]

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef int i_loop
        cdef double a1_pos[3],a2_pos[3],r12[3],d1[3],d2[3]
        cdef double lr12 = 0.0, lr32 = 0.0, r12_dot_r32 = 0.0, deriv
        coordinates = <vector3 *>input.coordinates.data

        for i_loop in range(3):
            a1_pos[i_loop] = coordinates[self.a1][i_loop]
            a2_pos[i_loop] = coordinates[self.a2][i_loop]
            r12[i_loop] = a2_pos[i_loop] - a1_pos[i_loop]
            if self.periodic:
                r12[i_loop] -= rint(r12[i_loop]/self.cellSize[i_loop])*self.cellSize[i_loop]
            lr12 += pow(r12[i_loop],2)
        lr12 = sqrt(lr12)
        dr = lr12 - self.r_eq

        #Add energy to MMTK, PIGS
        energy.energy_terms[self.index] = self.scale*self.D_r*(self.c1*pow((dr),2.0) - self.c2*pow(dr,3.0) + self.c3*pow(dr,4.0))
        
        #Add gradients
        if energy.gradients != NULL:
            deriv = -1.*self.D_r*(2.*self.c1*dr -3.*self.c2*pow(dr,2.0) + 4.*self.c3*pow(dr,3.0)) * self.scale #PIGS
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

            for i_loop in range(3):
                d1[i_loop] = r12[i_loop] / lr12
                d2[i_loop] = - d1[i_loop]

            for i_loop in range(3):
                gradients[self.a1][i_loop] += d1[i_loop]*deriv
                gradients[self.a2][i_loop] += d2[i_loop]*deriv

cdef class ElectrostaticTerm(EnergyTerm):


    cdef int atoms_per_molecule, periodic
    cdef double f,f2,electrostatic_energy, MM,HH,MH, cellSize[3]
    cdef atom_indices
    cdef scale

    def __init__(self, universe, atom_indices, frac, o_charge, e_energy, beg_end):
        cdef ArrayType ref_array
        self.atoms_per_molecule = len(atom_indices[0])
        self.atom_indices = atom_indices.flatten()
        self.f = frac
        self.f2 = 0.5 * (1. - frac)
        self.MM = 1.*o_charge**2
        self.HH = 1.*(o_charge * 0.5)**2
        self.MH = -1. * o_charge * (o_charge * 0.5)
        self.electrostatic_energy = e_energy
        self.periodic = universe.is_periodic
        self.scale = 1.
        if beg_end:
            self.scale = 0.5
        if self.periodic:
            self.cellSize[0] = universe.cellParameters()[0]
            self.cellSize[1] = universe.cellParameters()[1]
            self.cellSize[2] = universe.cellParameters()[2]
        EnergyTerm.__init__(self, universe,
                            "electrostatic", ("electrostatic",))
        self.eval_func = <void *>ElectrostaticTerm.evaluate

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef int i, j, k
        cdef double energy_sum = 0.0
        cdef double m1[3],m2[3]
        cdef double h1h3[3],h1h4[3],h1m2[3],h2h3[3],h2h4[3],h2m2[3],m1h3[3],m1h4[3],m1m2[3]
        cdef double h1h3r2,h1h4r2,h1m2r2,h2h3r2,h2h4r2,h2m2r2,m1h3r2,m1h4r2,m1m2r2
        cdef double h1h3r,h1h4r,h1m2r,h2h3r,h2h4r,h2m2r,m1h3r,m1h4r,m1m2r, deriv, dx,dy,dz
        coordinates = <vector3 *>input.coordinates.data

        #Loop over all molecules except the last on3
        for i in range(0,len(self.atom_indices)-1,self.atoms_per_molecule):
            #Calculate the virtual charge site of molecule 1
            for k in range(3):
                m1[k] = self.f*coordinates[self.atom_indices[i+2]][k] + self.f2*coordinates[self.atom_indices[i]][k] + self.f2*coordinates[self.atom_indices[i+1]][k]
            #Loop over all molecules from i til the end
            for j in range(i+self.atoms_per_molecule,len(self.atom_indices),self.atoms_per_molecule):
                h1h3r2 = 0.
                h1h4r2 = 0.
                h1m2r2 = 0.
                h2h3r2 = 0.
                h2h4r2 = 0.
                h2m2r2 = 0.
                m1h3r2 = 0.
                m1h4r2 = 0.
                m1m2r2 = 0.
                #Calculate the virtual charge site of molecule 2
                for k in range(3):
                    m2[k] = self.f*coordinates[self.atom_indices[j+2]][k] + self.f2*coordinates[self.atom_indices[j]][k] + self.f2*coordinates[self.atom_indices[j+1]][k]
                    
                    h1h3[k] = coordinates[self.atom_indices[i]][k] - coordinates[self.atom_indices[j]][k]
                    h1h4[k] = coordinates[self.atom_indices[i]][k] - coordinates[self.atom_indices[j+1]][k]
                    h1m2[k] = coordinates[self.atom_indices[i]][k] - m2[k]
                    h2h3[k] = coordinates[self.atom_indices[i+1]][k] - coordinates[self.atom_indices[j]][k]
                    h2h4[k] = coordinates[self.atom_indices[i+1]][k] - coordinates[self.atom_indices[j+1]][k]
                    h2m2[k] = coordinates[self.atom_indices[i+1]][k] - m2[k]
                    m1h3[k] = m1[k] - coordinates[self.atom_indices[j]][k]
                    m1h4[k] = m1[k] - coordinates[self.atom_indices[j+1]][k]
                    m1m2[k] = m1[k] - m2[k]
                    
                    #Get nearest image if required
                    if self.periodic:
                        h1h3[k] -= rint(h1h3[k]/self.cellSize[k])*self.cellSize[k]
                        h1h4[k] -= rint(h1h4[k]/self.cellSize[k])*self.cellSize[k]
                        h1m2[k] -= rint(h1m2[k]/self.cellSize[k])*self.cellSize[k]
                        h2h3[k] -= rint(h2h3[k]/self.cellSize[k])*self.cellSize[k]
                        h2h4[k] -= rint(h2h4[k]/self.cellSize[k])*self.cellSize[k]
                        h2m2[k] -= rint(h2m2[k]/self.cellSize[k])*self.cellSize[k]
                        m1h3[k] -= rint(m1h3[k]/self.cellSize[k])*self.cellSize[k]
                        m1h4[k] -= rint(m1h4[k]/self.cellSize[k])*self.cellSize[k]
                        m1m2[k] -= rint(m1m2[k]/self.cellSize[k])*self.cellSize[k]
                        
                    h1h3r2 += pow(h1h3[k],2)
                    h1h4r2 += pow(h1h4[k],2)
                    h1m2r2 += pow(h1m2[k],2)
                    h2h3r2 += pow(h2h3[k],2)
                    h2h4r2 += pow(h2h4[k],2)
                    h2m2r2 += pow(h2m2[k],2)
                    m1h3r2 += pow(m1h3[k],2)
                    m1h4r2 += pow(m1h4[k],2)
                    m1m2r2 += pow(m1m2[k],2)

                h1h3r = pow(h1h3r2,0.5)
                h1h4r = pow(h1h4r2,0.5)
                h1m2r = pow(h1m2r2,0.5)
                h2h3r = pow(h2h3r2,0.5)
                h2h4r = pow(h2h4r2,0.5)
                h2m2r = pow(h2m2r2,0.5)
                m1h3r = pow(m1h3r2,0.5)
                m1h4r = pow(m1h4r2,0.5)
                m1m2r = pow(m1m2r2,0.5)

                #These energies will all be passed back together at the very end
                energy_sum += self.electrostatic_energy * self.HH / h1h3r
                energy_sum += self.electrostatic_energy * self.HH / h1h4r
                energy_sum += self.electrostatic_energy * self.MH / h1m2r
                energy_sum += self.electrostatic_energy * self.HH / h2h3r
                energy_sum += self.electrostatic_energy * self.HH / h2h4r
                energy_sum += self.electrostatic_energy * self.MH / h2m2r
                energy_sum += self.electrostatic_energy * self.MH / m1h3r
                energy_sum += self.electrostatic_energy * self.MH / m1h4r
                energy_sum += self.electrostatic_energy * self.MM / m1m2r

                if energy.gradients != NULL:
                    gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

                    #H1 H3 interaction
                    deriv = self.electrostatic_energy * self.HH / (h1h3r * h1h3r2) * self.scale
                    dx = deriv * h1h3[0]
                    dy = deriv * h1h3[1]
                    dz = deriv * h1h3[2]
                    gradients[self.atom_indices[i]][0] -= dx
                    gradients[self.atom_indices[j]][0] += dx
                    gradients[self.atom_indices[i]][1] -= dy
                    gradients[self.atom_indices[j]][1] += dy
                    gradients[self.atom_indices[i]][2] -= dz
                    gradients[self.atom_indices[j]][2] += dz

                    #H1 H4 interaction
                    deriv = self.electrostatic_energy * self.HH / (h1h4r * h1h4r2) * self.scale
                    dx = deriv * h1h4[0]
                    dy = deriv * h1h4[1]
                    dz = deriv * h1h4[2]
                    gradients[self.atom_indices[i]][0]   -= dx
                    gradients[self.atom_indices[j+1]][0] += dx
                    gradients[self.atom_indices[i]][1]   -= dy
                    gradients[self.atom_indices[j+1]][1] += dy
                    gradients[self.atom_indices[i]][2]   -= dz
                    gradients[self.atom_indices[j+1]][2] += dz

                    #H2 H3 interaction
                    deriv = self.electrostatic_energy * self.HH / (h2h3r * h2h3r2) * self.scale
                    dx = deriv * h2h3[0]
                    dy = deriv * h2h3[1]
                    dz = deriv * h2h3[2]
                    gradients[self.atom_indices[i+1]][0] -= dx
                    gradients[self.atom_indices[j]][0]   += dx
                    gradients[self.atom_indices[i+1]][1] -= dy
                    gradients[self.atom_indices[j]][1]   += dy
                    gradients[self.atom_indices[i+1]][2] -= dz
                    gradients[self.atom_indices[j]][2]   += dz

                    #H2 H4 interaction
                    deriv = self.electrostatic_energy * self.HH / (h2h4r * h2h4r2) * self.scale
                    dx = deriv * h2h4[0]
                    dy = deriv * h2h4[1]
                    dz = deriv * h2h4[2]
                    gradients[self.atom_indices[i+1]][0] -= dx
                    gradients[self.atom_indices[j+1]][0] += dx
                    gradients[self.atom_indices[i+1]][1] -= dy
                    gradients[self.atom_indices[j+1]][1] += dy
                    gradients[self.atom_indices[i+1]][2] -= dz
                    gradients[self.atom_indices[j+1]][2] += dz

                    #H1 M2 interaction, changes xyz of H1, H3,H4,O2
                    deriv = self.electrostatic_energy * self.MH / (h1m2r * h1m2r2) * self.scale
                    dx = deriv * h1m2[0]
                    dy = deriv * h1m2[1]
                    dz = deriv * h1m2[2]
                    gradients[self.atom_indices[i]][0] -= dx                #H1 x
                    gradients[self.atom_indices[i]][1] -= dy                #H1 y
                    gradients[self.atom_indices[i]][2] -= dz                #H1 z
                    gradients[self.atom_indices[j]][0] += dx*self.f2        #H3 x
                    gradients[self.atom_indices[j]][1] += dy*self.f2        #H3 y
                    gradients[self.atom_indices[j]][2] += dz*self.f2        #H3 z
                    gradients[self.atom_indices[j+1]][0] += dx*self.f2      #H4 x
                    gradients[self.atom_indices[j+1]][1] += dy*self.f2      #H4 y
                    gradients[self.atom_indices[j+1]][2] += dz*self.f2      #H4 z
                    gradients[self.atom_indices[j+2]][0] += dx*self.f       #O2 x
                    gradients[self.atom_indices[j+2]][1] += dy*self.f       #O2 y
                    gradients[self.atom_indices[j+2]][2] += dz*self.f       #O2 z

                    #H2 M2 interaction, changes xyz of H2, H3,H4,O2
                    deriv = self.electrostatic_energy * self.MH / (h2m2r * h2m2r2) * self.scale
                    dx = deriv * h2m2[0]
                    dy = deriv * h2m2[1]
                    dz = deriv * h2m2[2]
                    gradients[self.atom_indices[i+1]][0] -= dx              #H2 x
                    gradients[self.atom_indices[i+1]][1] -= dy              #H2 y
                    gradients[self.atom_indices[i+1]][2] -= dz              #H2 z
                    gradients[self.atom_indices[j]][0] += dx*self.f2        #H3 x
                    gradients[self.atom_indices[j]][1] += dy*self.f2        #H3 y
                    gradients[self.atom_indices[j]][2] += dz*self.f2        #H3 z
                    gradients[self.atom_indices[j+1]][0] += dx*self.f2      #H4 x
                    gradients[self.atom_indices[j+1]][1] += dy*self.f2      #H4 y
                    gradients[self.atom_indices[j+1]][2] += dz*self.f2      #H4 z
                    gradients[self.atom_indices[j+2]][0] += dx*self.f       #O2 x
                    gradients[self.atom_indices[j+2]][1] += dy*self.f       #O2 y
                    gradients[self.atom_indices[j+2]][2] += dz*self.f       #O2 z

                    #M1 H3 interaction, changes xyz of H3, H1,H2,O1
                    deriv = self.electrostatic_energy * self.MH / (m1h3r * m1h3r2) * self.scale
                    dx = deriv * m1h3[0]
                    dy = deriv * m1h3[1]
                    dz = deriv * m1h3[2]
                    gradients[self.atom_indices[i]][0] -= dx*self.f2        #H1 x
                    gradients[self.atom_indices[i]][1] -= dy*self.f2        #H1 y
                    gradients[self.atom_indices[i]][2] -= dz*self.f2        #H1 z
                    gradients[self.atom_indices[i+1]][0] -= dx*self.f2      #H2 x
                    gradients[self.atom_indices[i+1]][1] -= dy*self.f2      #H2 y
                    gradients[self.atom_indices[i+1]][2] -= dz*self.f2      #H2 z
                    gradients[self.atom_indices[i+2]][0] -= dx*self.f       #O1 x
                    gradients[self.atom_indices[i+2]][1] -= dy*self.f       #O1 y
                    gradients[self.atom_indices[i+2]][2] -= dz*self.f       #O1 z
                    gradients[self.atom_indices[j]][0] += dx                #H3 x
                    gradients[self.atom_indices[j]][1] += dy                #H3 y
                    gradients[self.atom_indices[j]][2] += dz                #H3 z

                    #M1 H4 interaction, changes xyz of H4, H1,H2,O1
                    deriv = self.electrostatic_energy * self.MH / (m1h4r * m1h4r2) * self.scale
                    dx = deriv * m1h4[0]
                    dy = deriv * m1h4[1]
                    dz = deriv * m1h4[2]
                    gradients[self.atom_indices[i]][0] -= dx*self.f2        #H1 x
                    gradients[self.atom_indices[i]][1] -= dy*self.f2        #H1 y
                    gradients[self.atom_indices[i]][2] -= dz*self.f2        #H1 z
                    gradients[self.atom_indices[i+1]][0] -= dx*self.f2      #H2 x
                    gradients[self.atom_indices[i+1]][1] -= dy*self.f2      #H2 y
                    gradients[self.atom_indices[i+1]][2] -= dz*self.f2      #H2 z
                    gradients[self.atom_indices[i+2]][0] -= dx*self.f       #O1 x
                    gradients[self.atom_indices[i+2]][1] -= dy*self.f       #O1 y
                    gradients[self.atom_indices[i+2]][2] -= dz*self.f       #O1 z
                    gradients[self.atom_indices[j+1]][0] += dx              #H4 x
                    gradients[self.atom_indices[j+1]][1] += dy              #H4 y
                    gradients[self.atom_indices[j+1]][2] += dz              #H4 z

                    #M1 M2 interaction, changes xyz of H4, H1,H2,O1
                    deriv = self.electrostatic_energy * self.MM / (m1m2r * m1m2r2) * self.scale
                    dx = deriv * m1m2[0]
                    dy = deriv * m1m2[1]
                    dz = deriv * m1m2[2]
                    gradients[self.atom_indices[i]][0] -= dx*self.f2        #H1 x
                    gradients[self.atom_indices[i]][1] -= dy*self.f2        #H1 y
                    gradients[self.atom_indices[i]][2] -= dz*self.f2        #H1 z
                    gradients[self.atom_indices[i+1]][0] -= dx*self.f2      #H2 x
                    gradients[self.atom_indices[i+1]][1] -= dy*self.f2      #H2 y
                    gradients[self.atom_indices[i+1]][2] -= dz*self.f2      #H2 z
                    gradients[self.atom_indices[i+2]][0] -= dx*self.f       #O1 x
                    gradients[self.atom_indices[i+2]][1] -= dy*self.f       #O1 y
                    gradients[self.atom_indices[i+2]][2] -= dz*self.f       #O1 z
                    gradients[self.atom_indices[j]][0] += dx*self.f2        #H3 x
                    gradients[self.atom_indices[j]][1] += dy*self.f2        #H3 y
                    gradients[self.atom_indices[j]][2] += dz*self.f2        #H3 z
                    gradients[self.atom_indices[j+1]][0] += dx*self.f2      #H4 x
                    gradients[self.atom_indices[j+1]][1] += dy*self.f2      #H4 y
                    gradients[self.atom_indices[j+1]][2] += dz*self.f2      #H4 z
                    gradients[self.atom_indices[j+2]][0] += dx*self.f       #O2 x
                    gradients[self.atom_indices[j+2]][1] += dy*self.f       #O2 y
                    gradients[self.atom_indices[j+2]][2] += dz*self.f       #O2 z
        
        #Send final energy back to MMTK
        energy.energy_terms[self.index] = self.scale*energy_sum


cdef class LennardJonesTerm(EnergyTerm):

    cdef int a1,a2, periodic
    cdef double A, B,cellSize[3]
    cdef atom_indices
    cdef scale

    def __init__(self, universe,atom_indices, epsilon, sigma, beg_end):
        cdef ArrayType ref_array
        cdef int i
        EnergyTerm.__init__(self, universe,
                            "Lennard_Jones", ("Lennard_Jones",))
        self.atom_indices = atom_indices
        eps = epsilon
        sig = sigma
        self.A = 4. * epsilon * pow(sigma,12)
        self.B = 4. * epsilon * pow(sigma,6)
        self.periodic = universe.is_periodic
        self.scale = 1.
        if beg_end:
            self.scale = 0.5
        if self.periodic:
            self.cellSize[0] = universe.cellParameters()[0]
            self.cellSize[1] = universe.cellParameters()[1]
            self.cellSize[2] = universe.cellParameters()[2]
        self.eval_func = <void *>LennardJonesTerm.evaluate

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef int i,j,atoms
        cdef double a1_pos[3],a2_pos[3],dx,dy,dz,r2,r2i,r6i,r12i,energy_sum=0.0
        cdef double lr12 = 0.0, lr32 = 0.0, r12_dot_r32 = 0.0, deriv, delta_x,delta_y,delta_z
        atoms = len(self.atom_indices)
        coordinates = <vector3 *>input.coordinates.data
       
        for i in range(atoms):
            for j in range(i+1,atoms):
                dx = coordinates[self.atom_indices[i]][0]-coordinates[self.atom_indices[j]][0]
                dy = coordinates[self.atom_indices[i]][1]-coordinates[self.atom_indices[j]][1]
                dz = coordinates[self.atom_indices[i]][2]-coordinates[self.atom_indices[j]][2]
                if self.periodic:
                    dx -= rint(dx/self.cellSize[0])*self.cellSize[0]
                    dy -= rint(dy/self.cellSize[1])*self.cellSize[1]
                    dz -= rint(dz/self.cellSize[2])*self.cellSize[2]
                r2 = dx*dx + dy*dy + dz*dz
                r2i = 1./r2
                r6i = r2i*r2i*r2i
                r12i = r6i*r6i

                #Total energy will be passed back at the very end
                energy_sum += (self.A * r12i) - (self.B * r6i)

                if energy.gradients != NULL:
                    deriv = (12*self.A*r12i - 6.*self.B*r6i) * r2i * self.scale
                    gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

                    delta_x = deriv * dx
                    delta_y = deriv * dy
                    delta_z = deriv * dz
                    gradients[self.atom_indices[i]][0] -= delta_x
                    gradients[self.atom_indices[j]][0] += delta_x
                    gradients[self.atom_indices[i]][1] -= delta_y
                    gradients[self.atom_indices[j]][1] += delta_y
                    gradients[self.atom_indices[i]][2] -= delta_z
                    gradients[self.atom_indices[j]][2] += delta_z

        #Return total energy
        energy.energy_terms[self.index] = self.scale*energy_sum
