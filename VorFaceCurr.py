## Highlight Voronoi Faces
#This modifier function lets you select and highlight specific faces of Voronoi polyhedra
#based on their Voronoi order.
from ovito.data import *
from ovito.modifiers import AssignColorModifier, DeleteSelectedModifier
from traits.api import *
from ovito.pipeline import ModifierInterface
from ovito.traits import ColorTrait
class HighlightVoroFaces(ModifierInterface):
    voronoi_order = Int(4, label = "Voronoi Order")
    voronoi_surfaceArea = Float(1, label = "Face Surface Area")
    color = ColorTrait(color = "yellow", label = "Face Color")
    color_bonds = Bool(label="Apply color to bonds")
    def modify(self, data, **kwargs):
        faces = data.surfaces['voronoi-polyhedra_'].faces_
        particles = data.particles
        selection = (faces['Voronoi Order'] == 4)
        faces.create_property('Highlighted', dtype=int, data=selection)
        if 'voronoi-polyhedra' not in data.surfaces:
            raise RuntimeError("Please first generate Voronoi polyhedra and bonds by using the Voronoi Analyis Modifier.")
        faces = data.surfaces['voronoi-polyhedra_'].faces_
        #Highlight faces with chosen Voronoi order
        sel = faces.create_property("Selection")
        sel[faces['Voronoi Order'] == 19293] = 1 #change to self.voronoi_order to work
        assign_color_mod = AssignColorModifier(operate_on = 'surface_faces:voronoi-polyhedra/faces', color=self.color, keep_selection = False)
        data.apply(assign_color_mod)
        if self.color_bonds:
            assign_color_mod.operate_on = 'bonds'
            data.apply(assign_color_mod)
            #Show only Voronoi bonds associated with highlighted faces
        bond_indices_to_delete = faces['Bond Index'][faces['Voronoi Order'] != self.voronoi_order]
        bond_del = faces['Bond Index'][faces['Area'] >= self.voronoi_surfaceArea]
        bond_sel = data.particles_.bonds_.create_property("Selection")
        bond_sel[bond_indices_to_delete] = 1
        bond_sel[bond_del] = 1
        data.apply(DeleteSelectedModifier(operate_on ={'bonds'}))
        print("This is a random assortment of face area for reference")
        print(faces['Area'][0:10])
        pdel = particles['Charge'] != 2
        psel = data.particles.create_property("Selection")
        psel[pdel] = 1
        data.apply(DeleteSelectedModifier(operate_on ={'particles'}))
        print(faces['Area'] >= 0.1)
        
        

