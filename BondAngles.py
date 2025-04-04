### Calculate Bond Angles of Particle
# This modifier outputs the angles between all pairwise combinations of bonds at one particle.

from ovito.data import *
import numpy as np 
from traits.api import Union, Int, Enum, Bool
from ovito.pipeline import ModifierInterface
from itertools import combinations

class CalculateBondAngles(ModifierInterface):
    
    mode = Enum("Particle Indices", "Particle Identifiers", label = "Operate on")
    particle_index = Union(None, Int, label="Particle")
    bond_output = Bool(label = "Output bonds", default = False)
    bond_mode = Enum("Bond Indices", "Bond Identifiers", label = "Output")

    def calculate_bond_angles(self, bond_vectors):
        b1 = bond_vectors[:,0]/np.linalg.norm(bond_vectors[:,0], axis = 1)[:, None]
        b2 = bond_vectors[:,1]/np.linalg.norm(bond_vectors[:,1], axis = 1)[:, None]
        return np.degrees(np.arccos(np.sum(b1*b2, axis = 1)))
    
    def calculate_bond_vector_combinations(self, data, particle):
        # Create bonds enumerator object.
        bonds_enum = BondsEnumerator(data.particles.bonds)
        # List of bond indices of all bonds at the current particle
        bonds_of_particle = [bond_index for bond_index in bonds_enum.bonds_of_particle(particle)]
        if len(bonds_of_particle) < 2:
            raise RuntimeError("Not enough bonds found to compute angles.")
        # All possible bond pairs in bonds_of_particle list
        bond_pairs = list(combinations(bonds_of_particle, 2))
        # Look up corresponding bond vectors
        topo = data.particles.bonds.topology[bond_pairs] 
        # Flip bond vector if current particle is not index 0 in topology pairs
        idx = np.where(topo[:,:,0] != particle)
        vectors = data.particles.bonds["Bond vectors"][bond_pairs]
        vectors[idx[0], idx[1]] *= -1
        # Get particle index triplets from topo [B,A][B,C]
        topo[idx[0], idx[1]] =  np.flip(topo[idx[0], idx[1]])
        triplets = np.column_stack((topo[:, 0, 1],topo[:, 1, :]))
        return bond_pairs, triplets, vectors
       
    def modify(self, data, **kwargs):
        if self.particle_index == None:
            return
        if self.mode == "Particle Identifiers" and "Particle Identifier" not in data.particles:
            raise RuntimeError("No Particle Identifiers in DataCollection. Deactivate Option <Output Particle Identifiers>.")
        if data.particles.bonds == None:
            raise RuntimeError("No Bonds in DataCollection. Please first generate bonds.")
        if self.bond_output and self.bond_mode == "Bond Identifiers" and "Bond Identifier" not in data.particles.bonds:
            raise RuntimeError("No Bond Identifiers in DataCollection. Switch to Option <Output Bond Indices>.")

        # Look up bond indices of all bonds connected to current particle
        if self.mode == 'Particle Identifiers':
            particle = np.where(data.particles.identifiers == self.particle_index)[0][0]
        else:
            particle = self.particle_index 

        # Calculate bond vectors for global topology array    
        positions = data.particles.positions
        topology = data.particles.bonds.topology
        bond_vectors = positions[topology[:,1]] - positions[topology[:,0]]
        bond_vectors += np.dot(data.cell[:3,:3], data.particles.bonds.pbc_vectors.T).T
        data.particles_.bonds_.create_property("Bond vectors", data = bond_vectors)    
            # Get all possible combinations of bond pairs at one particle, the corresponding 
        faces = data.surfaces['voronoi-polyhedra_'].faces_
        bond_VI = faces['Bond Index'][faces['Voronoi Order'] != 4]
        bond_del = faces['Bond Index'][faces['Area'] >= 0.3]
        bond_sel = data.particles_.bonds_.create_property("Selection")
        bond_sel[bond_VI] = 1
        bond_sel[bond_del] = 1
        Bonds = []
        Angles = []
        Bond_Pairs = []
        Particle_Triplets = []
                
        for j in range(data.particles.count):
            yield j/data.particles.count
            # particle index triplets and bond vectors
            bond_pairs, triplets, v_b = self.calculate_bond_vector_combinations(data, j)
            
            # Calculcate angles between all pairs of bond vectors
            if v_b is None:
                continue
            
            Angles+=list(self.calculate_bond_angles(v_b))
            Bond_Pairs+=bond_pairs
            Particle_Triplets+=list(triplets)
    
        # Store results as OVITO Data Table#
        for a in Bond_Pairs:
            Bonds.append([bond_sel[a[0]], bond_sel[a[1]]])
            
        Angles = np.reshape(Angles, -1)
        Bond_Pairs = np.reshape(np.array(Bond_Pairs), (-1, 2))
        Triplets = np.reshape(np.array(Particle_Triplets), (-1, 3))
        table = data.tables.create(
            identifier="bond-angles",
            title=f"Bond",
            plot_mode=DataTable.PlotMode.NoPlot)
        table.y = table.create_property('Angle', data=Angles)
        table.y = table.create_property('Particle Triplet', data=Triplets, components=['A', 'B', 'C'])
        table.y = table.create_property('Bond Pair', data=Bond_Pairs, components=['Bond1', 'Bond2'])
        table.y = table.create_property('Bond VI', data=Bonds)
        
      
           

       
      

            
