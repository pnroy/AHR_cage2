#!/usr/bin/env python
# -*- coding: utf-8 -*-

# q-TIP4P water potential
# Kevin Bishop

from MMTK.ForceFields.ForceField import ForceField
from MMTK_qTIP4p import HarmonicAngleTerm, HarmonicBondTerm, QuarticBondTerm, ElectrostaticTerm, LennardJonesTerm
from MMTK.Units import electrostatic_energy
import numpy as np

class HarmonicAngleForceField(ForceField):

    """
    Harmonic angle bend potential between 3 atoms
    """

    def __init__(self, atoms, theta_eq, k_theta):
        self.atom_indices = map(lambda x: self.getAtomParameterIndices([x])[0],atoms)
        self.theta_eq = theta_eq
        self.k_theta = k_theta
        self.arguments = (self.atom_indices,theta_eq,k_theta)
        ForceField.__init__(self, 'harmonic_angle')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")

        f, offsets = self.beadOffsetsAndFactor([self.atom_indices[0],self.atom_indices[1],self.atom_indices[2]], global_data)

        highest_offset = offsets[-1][-1]

        return [HarmonicAngleTerm(universe,
                                       self.atom_indices + o,
                                       self.theta_eq,
                                       self.k_theta,
                                       o[0]==0 or o[0]==highest_offset)  #True if first/last bead (Required for PIGS)
                for o in offsets]

class HarmonicBondForceField(ForceField):

    """
    Harmonic bond stretch potential between 2 atoms
    """

    def __init__(self, atoms, r_eq, k_r):
        self.atom_indices = map(lambda x: self.getAtomParameterIndices([x])[0],atoms)
        self.r_eq = r_eq
        self.k_r = k_r
        self.arguments = (self.atom_indices,r_eq,k_r)
        ForceField.__init__(self, 'harmonic_bond')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")

        f, offsets = self.beadOffsetsAndFactor([self.atom_indices[0],self.atom_indices[1]], global_data)

        highest_offset = offsets[-1][-1]

        return [HarmonicBondTerm(universe,
                                       self.atom_indices + o,
                                       self.r_eq,
                                       self.k_r,
                                       o[0]==0 or o[0]==highest_offset)
                for o in offsets]

class QuarticBondForceField(ForceField):

    """
    Quartic bond stretch potential between 2 atoms
    """

    def __init__(self, atoms, D_r, c1, c2, c3, r_eq):
        self.atom_indices = map(lambda x: self.getAtomParameterIndices([x])[0],atoms)
        self.D_r = D_r
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.r_eq = r_eq
        self.arguments = (self.atom_indices,D_r,c1,c2,c3,r_eq)
        ForceField.__init__(self, 'quartic_bond')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):

        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        f, offsets = self.beadOffsetsAndFactor([self.atom_indices[0],self.atom_indices[1]], global_data)

        highest_offset = offsets[-1][-1]

        return [QuarticBondTerm(universe,
                                       self.atom_indices + o,
                                       self.D_r,
                                       self.c1,
                                       self.c2,
                                       self.c3,
                                       self.r_eq,
                                       o[0]==0 or o[0]==highest_offset)
                for o in offsets]

class ElectrostaticForceField(ForceField):

    """
    Electrostatic interactions between all charged particles in universe
    WILL ONLY WORK WITH WATERS! 3 AND 4 POINT MODELS SPECIFICALLY!
    """

    def __init__(self, universe,fraction,o_charge):
        #Take the list of objects (molecules) and create a new list where each element in the
        #list is another list containing the indices of all the atoms in the molecule
        waters = universe.objectList()
        water_atoms = map(lambda x: map(lambda y: self.getAtomParameterIndices([y])[0],x.atomList()), waters)
        self.atom_indices = water_atoms
        self.fraction = fraction
        self.o_charge = o_charge
        self.arguments = (self.atom_indices,self.fraction,self.o_charge)
        ForceField.__init__(self, 'electrostatics')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):

        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")

        f, offsets = self.beadOffsetsAndFactor(map(lambda x: x[0],self.atom_indices), global_data)
        
        highest_offset = offsets[-1][-1]

        return [ElectrostaticTerm(universe,
                                       np.array(map(lambda x: self.atom_indices[x]+o[x],range(len(o)))),
                                       self.fraction,
                                       self.o_charge,
                                       electrostatic_energy,
                                       o[0]==0 or o[0]==highest_offset)
                for o in offsets]

class LennardJonesForceField(ForceField):

    """
    Lennard-Jones potential between all of the Oxygens in the water molecules
    """

    def __init__(self, universe, epsilon, sigma):
        #Get atom indices for oxygens
        waters = universe.objectList()
        oxygens = map(lambda x: x.atomList()[2],waters)     #Oxygen is 2 in atom list
        atom_indices = map(lambda x: self.getAtomParameterIndices([x])[0],oxygens)
        self.atom_indices = atom_indices
        self.epsilon = epsilon
        self.sigma = sigma

        self.arguments = (self.atom_indices,self.epsilon,self.sigma)
        ForceField.__init__(self, 'lennard_jones')

    def ready(self, global_data):
        return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        
        f, offsets = self.beadOffsetsAndFactor(self.atom_indices, global_data)

        highest_offset = offsets[-1][-1]

        return [LennardJonesTerm(universe,
                                       self.atom_indices + o,
                                       self.epsilon,
                                       self.sigma,
                                       o[0]==0 or o[0]==highest_offset)
                for o in offsets]
                