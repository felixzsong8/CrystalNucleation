##Code for Ovito Data
#from ovito.io import import_file
from tqdm import trange
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#
#
import freud
import gsd.hoomd
import matplotlib.pyplot as plt
 

#Read through the simulation data
#Look at each particle to decide if its either crystal bulk, crystal interface, liquid interface or liquid bulk
pnum = 2048
fnum = 200
infile = "/Users/felixsong/Desktop/Ovito/mc55_2048_150.gsd" # your file path 
traj = gsd.hoomd.open(infile, 'rb')
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
interfaceP = []

for framedata in traj:
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    qlnc = ql.num_connections
    nnlist = []
    for i,j in ql.nlist:
        nnlist.append([i,j])
    nlistdf = pd.DataFrame(nnlist,columns=['i','j'])
    s = []
    pair = [0,0]
    count = 0
    for b in range(len(nlistdf.j)):
        if(qlnc[nlistdf.j[b]] < 8): #Saying this neighbor particle is liquid
            pair[0] += 1
        else:
            pair[1] += 1 #Saying this neighbor particle is crystal
        #reset for each particle
        count += 1
        if(count == 12):
            if(qlnc[nlistdf.i[b]] < 8):
                if(pair[1] != 0):
                    s.append(2) #liquid interface
                else:
                    s.append(3) #liquid bulk
            else:
                if(pair[0] != 0):
                    s.append(1) #crystal interface
                else:
                    s.append(0) #crystal in bulk
            pair = [0,0] 
            count = 0
    
    interfaceP.append(s)
    

#%%

def check(c, nlist, array):
    for a in nlist:
        if(qlnc[a] >= 8 and qlci[a] != qlci[c]):
            return array[a]
    return qlci[c]
#%%
#Making clusterId as an array with cluster ids
#the liquid interface partilces are set a -2, liquid bulk at -1
#Single crystals are set at -5 and go down
#Clusters are set at 0 and go up

infile = "/Users/felixsong/Desktop/Ovito/sim/cluster.gsd" # your file path 
traj = gsd.hoomd.open(infile, 'rb')
clusterId = []
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
frame = 0
idxCap = int(2048/100)
for framedata in traj:
    array = []
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    qlnc = ql.num_connections
    qlci = ql.cluster_idx
    nnlist = []
    for i,j in ql.nlist:
        nnlist.append([i,j])
    nlistdf = pd.DataFrame(nnlist,columns=['i','j'])
    sCount = -5
    for c in range(pnum):
        if(framedata.particles.charge[c] == 3): #Liquid particles 
            array.append(-1)
        elif(framedata.particles.charge[c] == 2): #Liquid interface 
            array.append(-2)
        else: #check to see if we are looking at single particle, which should have no crystal neighbors but be crystal
            count = 12*c
            single = 0
            for d in range(12):
                if(framedata.particles.charge[nlistdf.j[count+d]] == 1 or framedata.particles.charge[nlistdf.j[count+d]] == 0 ):
                    single += 1
            if(single == 0):
                array.append(sCount)
                sCount -= 1
            else:
                array.append(qlci[c])
    #After indexing all of the partilces, check to make sure all crystal particles
    #with a clusterindex too large dont have any near neighbors. If they don't, treat as single
    for c in range(pnum):
        if(array[c] > idxCap and framedata.particles.charge[c] <= 1):
            value = check(c, nlistdf.j[c*12:c*12+11], array)
            if(value >= idxCap):
                array[c] = sCount 
                sCount -= 1
            else:
                array[c] = value
            print("corrected")
    clusterId.append(array)
    frame += 1
    print(frame)
    
    




#%%

#Making a gsd File based on previous gsd file, including new parameter, charge


setnum = 150
particlenum = 2048
#infile = '/Users/felixsong/Desktop/Ovito/mc55_2048_150.gsd'
outfile = '/Users/felixsong/Desktop/Ovito/sim/clusterID.gsd' #Change this to cluster.gsd and change charge property to interfaceP if you want to not have the id
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
traj = gsd.hoomd.open(infile, 'rb')
framenum = len(traj)
trajout = gsd.hoomd.open(outfile, 'wb') 
count = 0 
for framedata in (traj):
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    s = gsd.hoomd.Snapshot()
    s.particles.N = particlenum
    s.particles.position = framedata.particles.position
    s.configuration.box = [12.4927,12.4927,12.4927,0,0,0]
    s.particles.charge = clusterId[count]
    #s.particles.charge = clusterId[count*pnum:(count+1)*pnum]
    count += 1
#     s.particles.mass = tempdf['cluster_idx78pp']
    trajout.append(s)

#%%

#Making a dataframe from a gsd file

import freud
import gsd.hoomd
import matplotlib.pyplot as plt
#Making a dataframe file
# set up cluster size change threshold
# need to rewrite, make clustersize a property to all the particles
change_threshold = 5 # or start with 3
qla = freud.order.Steinhardt(l=6,average=True,wl=False)
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
cl_props = freud.cluster.ClusterProperties()
particlenum = pnum
framenum = fnum
growlist = []
# dataray = np.zeros((1045,200,2048,4))
datalist = []

infile = '/Users/felixsong/Desktop/Ovito/sim/cluster.gsd'
traj = gsd.hoomd.open(infile, 'rb')
   
for frame in traj:
    ql.compute(system=frame, neighbors={'num_neighbors': 12})
    qla.compute(system=frame, neighbors={'num_neighbors': 12})
    cl_props.compute(frame, ql.cluster_idx)
    size = ql.cluster_sizes[ql.cluster_idx]
    rg = cl_props.radii_of_gyration[ql.cluster_idx]
    aveql = qla.particle_order
    position = frame.particles.position
    s = frame.particles.charge
    
    for i in range(particlenum):
        datalist.append([size[i],rg[i],aveql[i],frame,i,s[i],position])


df1 = pd.DataFrame(datalist,columns=['size1','rg','aveql','frame','particle','ParticleState','Position']) 


#%%
    
    