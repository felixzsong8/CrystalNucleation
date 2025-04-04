#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:28:22 2024

@author: felixsong
"""

## Highlight Voronoi Faces
#This modifier function lets you select and highlight specific faces of Voronoi polyhedra
#based on their Voronoi order.
from ovito.data import *
from ovito.modifiers import AssignColorModifier, DeleteSelectedModifier
from traits.api import *
from ovito.pipeline import ModifierInterface
from ovito.traits import ColorTrait
import numpy as np
from traits.api import Union, Int, Enum, Bool
from itertools import combinations
from functools import reduce
pnum = 2048
partDist = [0] * pnum
class CalculateBondAngles(ModifierInterface):
    mode = Enum("Particle Indices", "Particle Identifiers", label = "Operate on")
    particle_index = Union(None, Int, label="Particle")
    bond_output = Bool(label = "Output bonds", default = False)
    bond_mode = Enum("Bond Indices", "Bond Identifiers", label = "Output")
    particle_index = 18
    def calculate_bond_angles(self, bond_vectors):
        b1 = bond_vectors[:,0]/np.linalg.norm(bond_vectors[:,0], axis = 1)[:, None]
        b2 = bond_vectors[:,1]/np.linalg.norm(bond_vectors[:,1], axis = 1)[:, None]
        return np.degrees(np.arccos(np.sum(b1*b2, axis = 1)))
    def modify(self, data, **kwargs):
        if self.particle_index == None:
            return
        if self.mode == "Particle Identifiers" and "Particle Identifier" not in data.particles:
            raise RuntimeError("No Particle Identifiers in DataCollection. Deactivate Option <Output Particle Identifiers>.")
        if data.particles.bonds == None:
            raise RuntimeError("No Bonds in DataCollection. Please first generate bonds.")
        if self.bond_output and self.bond_mode == "Bond Identifiers" and "Bond Identifier" not in data.particles.bonds:
            raise RuntimeError("No Bond Identifiers in DataCollection. Switch to Option <Output Bond Indices>.")
        #Compute bond vectors
        topology = data.particles.bonds.topology
        positions = data.particles.positions
        bond_vectors = positions[topology[:,1]] - positions[topology[:,0]]
        bond_vectors += np.dot(data.cell[:3,:3], data.particles.bonds.pbc_vectors.T).T
        faces = data.surfaces['voronoi-polyhedra_'].faces_    
        data.particles_.bonds_.create_property("Bond vectors", data = bond_vectors)
        # Create bonds enumerator object.
        bonds_enum = BondsEnumerator(data.particles.bonds)
        # Look up bond indices of all bonds connected to current particle
        bonds_of_particle = [bond_index for bond_index in bonds_enum.bonds_of_particle(self.particle_index)]
        # All possible bond pairs in bonds_of_particle list
        bond_pairs = list(combinations(bonds_of_particle, 2))
        if len(bond_pairs) > 1:
            angles = self.calculate_bond_angles(np.reshape(bond_vectors[bond_pairs], (len(bond_pairs), 2, -1)))
            header = "Triplet A - B - C \t Angle"
            if self.bond_output:
                header += "\t Bond Pair B1 - B2"
            print(header)
            print(f"{'-'*(len(header)+16)}")
            b = len(str(data.particles.bonds.count))
            p = len(str(data.particles.count))
            faces = data.surfaces['voronoi-polyhedra_'].faces_
            bond_VI = faces['Bond Index'][faces['Voronoi Order'] != 4]
            bond_sel = data.particles_.bonds_.create_property("Selection")
            bond_sel[bond_VI] = 1
            Triplets = []
            Bonds = []
            for i in range(len(bond_pairs)):
                triplet = topology[[bond_pairs[i]]]
                B = self.particle_index
                A, C = triplet[triplet!=B]
                if self.mode == 'Particle Identifiers':
                    A,B,C = data.particles.identifiers[[A, self.particle_index, C]]
                Triplets.append([A,B,C])
                output = f"{A:{p}d} - {B} - {C:<{p}d} \t {angles[i]:7.3f}"
                if self.bond_output:
                    B1, B2 = bond_pairs[i][0], bond_pairs[i][1]
                    if self.bond_mode == 'Bond Identifiers':
                        B1, B2 = int(data.particles.bonds['Bond Identifier'][B1]),int(data.particles.bonds['Bond Identifier'][B2])
                    output+=f"\t {B1:{b}d} - {B2:<{b}d}"
                    Bonds.append([bond_sel[B1],bond_sel[B2]])        
                print(output)
            table = data.tables.create(
                identifier="bond-angles",
                title=f"Bond Angles of Particle {self.particle_index}",
                plot_mode=DataTable.PlotMode.NoPlot)
            table.create_property('Angle', data=angles)
            table.create_property('Particle Triplet', data=Triplets, components=['A', 'B', 'C'])
            # table.create_property('Voronoi Order', data=bond_VI)
            table.create_property('Bond VI', data=Bonds)
            #Make a new property here for the voronoi order
            if self.bond_output:
                if self.bond_mode == 'Bond Identifiers':
                    bond_pairs = data.particles.bonds['Bond Identifier'][bond_pairs]
                table.create_property('Bond Pair', data=bond_pairs, components=['Bond1', 'Bond2'])
            for a in range(len(Bonds)):
                if(Bonds[a][0] == 0 and Bonds[a][1] == 0):
                    partDist[self.particle_index] += 1
            print(partDist[10:20])
        else:
            print("Not enough bonds found to compute angles.")
