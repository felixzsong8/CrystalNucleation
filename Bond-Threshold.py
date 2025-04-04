## Highlight Voronoi Faces
#This modifier function lets you select and highlight specific faces of Voronoi polyhedra
#based on their Voronoi order.
from ovito.data import *
from ovito.modifiers import AssignColorModifier, DeleteSelectedModifier
from traits.api import *
from ovito.pipeline import ModifierInterface
from ovito.traits import ColorTrait
class HighlightVoroFaces(ModifierInterface):
    voronoi_order = Int(6, label = "Voronoi Order")
    color = ColorTrait(color = "yellow", label = "Face Color")
    color_bonds = Bool(label="Apply color to bonds")
    def modify(self, data, **kwargs):
        faces = data.surfaces['voronoi-polyhedra_'].faces_
        selection = (faces['Voronoi Order'] == 4)
        faces.create_property('Highlighted', dtype=int, data=selection)
        if 'voronoi-polyhedra' not in data.surfaces:
            raise RuntimeError("Please first generate Voronoi polyhedra and bonds by using the Voronoi Analyis Modifier.")
        faces = data.surfaces['voronoi-polyhedra_'].faces_
        #Highlight faces with chosen Voronoi order
        sel = faces.create_property("Selection")
        sel[faces['Voronoi Order'] == self.voronoi_order] = 1
        assign_color_mod = AssignColorModifier(operate_on = 'surface_faces:voronoi-polyhedra/faces', color=self.color, keep_selection = False)
        data.apply(assign_color_mod)
        if self.color_bonds:
            assign_color_mod.operate_on = 'bonds'
            data.apply(assign_color_mod)
            #Show only Voronoi bonds associated with highlighted faces
        bond_indices_to_delete = faces['Bond Index'][faces['Voronoi Order'] != self.voronoi_order]
        bond_del = faces['Bond Index'][faces['Area'] >= 1]
        bond_sel = data.particles_.bonds_.create_property("Selection")
        bond_sel[bond_indices_to_delete] = 1
        bond_sel[bond_del] = 1
        data.apply(DeleteSelectedModifier(operate_on ={'bonds'}))
        print(bond_sel[0:20])
        print(faces['Area'][0:10])
        print(bond_del[0:10])
        print(bond_indices_to_delete[0:10])

