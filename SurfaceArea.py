##Highlight Voronoi Faces
#This modifier function lets you select and highlight specific faces of Voronoi polyhedra
#based on their Voronoi order.
from ovito.data import *
from ovito.modifiers import AssignColorModifier, DeleteSelectedModifier
import numpy
def modify(frame: int, data: DataCollection):
    areaArray = []
    for a in range(10):#This is for number of clusters, if you think there are more, increase
        areaArray.append(0)
        mesh = data.surfaces['voronoi-polyhedra_']
        ptypes = data.particles['Charge']
        sel = mesh.faces_.create_property("Selection", data = numpy.zeros(mesh.faces.count))
        
        A = ptypes[mesh.faces["Region"]]
        B = ptypes[mesh.faces["Adjacent Cell"]]
   
        sel[:] = numpy.logical_not(numpy.logical_and(A == a, B >= -9))
        
       # mesh.delete_faces( mesh.faces['Selection'] == 1) 
        value = numpy.sum(mesh.faces["Area"][mesh.faces['Selection'] != 1])
        areaArray[a] += value
    print("Cluster Area")
    print(areaArray)
    singleArray = []
    for a in range(8):#This is for clusters
        singleArray.append(0)
        mesh = data.surfaces['voronoi-polyhedra_']
        ptypes = data.particles['Charge']
        sel = mesh.faces_.create_property("Selection", data = numpy.zeros(mesh.faces.count))
        
        A = ptypes[mesh.faces["Region"]]
        B = ptypes[mesh.faces["Adjacent Cell"]]
   
        sel[:] = numpy.logical_not(numpy.logical_and(A == -1-a, B >= -9))
        
       # mesh.delete_faces( mesh.faces['Selection'] == 1) 
        value = numpy.sum(mesh.faces["Area"][mesh.faces['Selection'] != 1])
        singleArray[a] += value
    
    
    print("Single Crystal Area")
    print(singleArray)
    clusterId = data.particles['Charge']
    clusterArea = []
    for a in range(2048):
        if(clusterId[a] >= 0):
            clusterArea.append(areaArray[int(clusterId[a])])
        else:
            clusterArea.append(0)
    table = data.tables.create(
        identifier="SurfaceArea",
        title="Area",
        plot_mode=DataTable.PlotMode.NoPlot)
    table.y = table.create_property('ClusterID', data=clusterId)
    table.y = table.create_property('ClusterArea', data=clusterArea)
        
        
    
    
        
        
        

