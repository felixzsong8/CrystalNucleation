##Code for Ovito Data
#from ovito.io import import_file
from tqdm import trange
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

filename = "/Users/felixsong/Desktop/Ovito/mc55_2048_150"
clist = ['x','y','z','type','ix','iy','iz','Coordination','AtomicVolume','CavityRadius','VoronoiIndex1','VoronoiIndex2','VoronoiIndex3','VoronoiIndex4','VoronoiIndex5','VoronoiIndex6','VoronoiIndex7','MaxFaceOrder','id']
#clist = ['x','y','z','type', 'ix', 'iy', 'iz', 'Coordination', 'id', 'StructureType', 'ColorR', 'ColorG', 'ColorB']
#This is for melt687001
#clist = ['x','y','z','type','ix','iy','iz','Coordination','id','AtomicVolume','CavityRadius','VoronoiIndex1','VoronoiIndex2','VoronoiIndex3','VoronoiIndex4','VoronoiIndex5','VoronoiIndex6','MaxFaceOrder']
#This is for relax484
#clist = ['x','y','z','type','ix','iy','iz','Coordination','AtomicVolume','CavityRadius','VoronoiIndex1','VoronoiIndex2','VoronoiIndex3','VoronoiIndex4','VoronoiIndex5','VoronoiIndex6','MaxFaceOrder','id']
def readOvito(filename,columnnames,downcastlist=[]):
    line = open(filename,'r').read().split('\n')
    particlenum = int(line[3]) # number of atoms
    framenum = len(line)//(particlenum+9)
    flist = []
    for frame in trange(framenum):
        data = line[frame*(particlenum+9)+9:(frame+1)*(particlenum+9)]
        for i in range(particlenum):
            floats = [float(x) for x in data[i].split()] # extract data
            floats.append(frame) # append framenum
            flist.append(floats)

    df = pd.DataFrame(flist, columns = columnnames + ['frame'])
    for i in downcastlist:
        df[columnnames[i]] = pd.to_numeric(df[columnnames[i]],downcast='integer')
    return df, particlenum, framenum



df, pnum, fnum = readOvito(filename, clist)
df['volumefraction'] = np.pi/6 / df['AtomicVolume']
#print(df)

#df.plot(kind = 'scatter', x = 'frame', y = 'a14', color = 'red')
#plt.show()


#%% Making data to plot 2d graph atomic volume to index4 
#The row Length should be number of frames
#The volumn length should be the number of particles

volume = df.AtomicVolume
voronoiIndex4 = df.VoronoiIndex4
row_length = fnum
column_length = pnum
volumeList = [[0 for j in range(column_length)] for i in range(row_length)]
index4List = [[0 for j in range(column_length)] for i in range(row_length)]
count = 0





for b in range(row_length):
    for c in range(column_length):
        index = (b * column_length) + c
        volumeList[b][c] = volume[index]
        
for d in range(row_length):
    for e in range(column_length):
        index = (d * column_length) + e
        index4List[d][e] = voronoiIndex4[index]
    
    
    
#%%

#This is for a total violin plot where all data is on one plot
sns.violinplot(data=df, x=df.drop(count).VoronoiIndex4, y=df.drop(count).volumefraction).set(title="mc55_2048_150: VoronoiIndex4 against volumefraction")
#sns.violinplot(data=df, x=df.drop(count).VoronoiIndex4[0:200], y=df.drop(count).AtomicVolume[0:200]).set(title="mc55_size10_2: VoronoiIndex4 against AtomicVolume")
#plt.ylim(0.3,0.7)
#plt.xlim(-0.5,6.5)
f = 38
#sns.violinplot(data=df, x=df.drop(count).VoronoiIndex4[pnum*f:pnum*f + pnum - 1], y=df.drop(count).volumefraction[pnum*f:pnum*f + pnum - 1]).set(title="relax484: VoronoiIndex4 against volumefraction " + str(f))


#%%
#Getting rid of the VoronoiIndex4 == 9, might be skewing data
#Add a or statement or whatever elements you don't want to include
count = []
for f in range(len(df.frame)):
    if(df.VoronoiIndex4[f] == 9):
        #df.drop([f]) 
        count.append(f)
        



#%%
#To plot the data frame by frame, using violin plots


for f in range(fnum):
    sns.violinplot(data=df, x=df.drop(count).VoronoiIndex4[pnum*f:pnum*f + pnum - 1], y=df.drop(count).volumefraction[pnum*f:pnum*f + pnum - 1]).set(title="relax484: VoronoiIndex4 against volumefraction " + str(f))
    #sns.violinplot(data=df, x=df.drop(count).VoronoiIndex4[pnum*f:pnum*f + pnum - 1], y=df.drop(count).AtomicVolume[pnum*f:pnum*f + pnum - 1]).set(title="mc55_size10_2: VoronoiIndex4 against AtomicVolume " + str(f))
    plt.ylim(0.3,0.8)
    plt.xlim(-0.5,6.5)
    plt.savefig("violinplot-" + str(f) + ".pdf")    
    plt.clf()
    
    

#%%
#Plots number of particles with certain parameters
#Index4 should have the length of the number of frames
index4 = [0] * fnum

for a in range(len(df.frame)):
    #if(df.VoronoiIndex4[a] >= 2) and (df.Coordination[a] >= 14):
    if(df.VoronoiIndex4[a] >= 3): #and (df.Coordination[a]>=14):
        index4[df.frame[a]] += 1
        


plt.plot(index4)
plt.title('mc55_2048_150: VoronoiIndex4 >= 3')
plt.show()

#%%
#Scatter plot of index to Atomic Volume
plt.scatter(volumeList, index4List, s=1)
plt.show()

print(len(volumeList[0]))

#%%

#plot scatter of number of particles with certain parameter

for f in range(fnum):
    array = [0] * 100
    for a in range(pnum):
        array[int(df.Coordination[pnum*f + a])] += 1
    plt.plot(array)
    #plt.ylim(0,1000)
    plt.xlim(5,25)
    plt.title("Num-melt-Coord, frame " + str(count))
    plt.savefig("particle-dist-" + str(f) + ".png")    
    plt.clf()
    
        
#%%        

plt.plot(array)
plt.show()


#%% 
#Loading information for Q6

import freud
import gsd.hoomd
import matplotlib.pyplot as plt

infile = "/Users/felixsong/Desktop/Ovito/melt.gsd" # your file path 
ql = freud.order.Steinhardt(l=6,average=True,wl=False) # try average= False as well
traj = gsd.hoomd.open(infile, 'rb')


#%%
#Plotting Q6 against volumeFraction

framenum = len(traj)
count = 0
for framedata in traj:
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    #ql.particle_order  this will be a list with length equal to total number of particles
    plt.scatter(ql.particle_order, df.volumefraction[count*pnum:count*pnum + pnum], s = 1)
    plt.ylim(0.3,0.7)
    plt.xlim(0,0.6)
    plt.title("Q6-volDensity-melt687, frame " + str(count))
    plt.savefig("Q6-volDensity-melt687-" + str(count) + ".png")    
    plt.clf()
    count += 1


#%%
#Plotting 3D plot between VolFrac, Q6, and VI4
from mpl_toolkits import mplot3d
fig = plt.figure()
framenum = len(traj)
count = 0
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
for framedata in traj:
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    ax = plt.axes(projection = '3d')
    xdata = df.volumefraction[count*pnum:count*pnum + pnum]
    zdata = df.VoronoiIndex4[count*pnum:count*pnum + pnum]
    ax.scatter3D(xdata, ql.num_connections, zdata, c = zdata, s=3)
    plt.xlim(0,0.8)
    plt.ylim(0,13)
    ax.set_zlim(0,10)
    ax.set(xlabel='VolumeFraction', ylabel='Q6-num-connections',zlabel='VoronoiIndex4')
    ax.set_title("3Dplot,VolFrac,Q6,VI4-melt-" + str(count))
    plt.savefig("Particle_dist-" + str(count) + ".png")    
    plt.clf()
    count += 1


#ql.particle_order
#ql.num_connections

#%%
#Function to create 

def num_coord(xdata,ydata,count,c):
    dict = {}
    array = [0] * len(xdata)

    for i in range(len(xdata)):
        if (xdata[i],ydata[count*pnum+i]) in dict:
            dict[xdata[i],ydata[count*pnum+i]] += 1
        else:
            dict[xdata[i],ydata[count*pnum+i]] = 1

    for j in range(len(xdata)):
        array[j] = c*dict[xdata[j],ydata[count*pnum+j]]
    return array

#%%

#2D plot
from collections import Counter



framenum = len(traj)
count = 0
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
for framedata in traj:
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    ydata = df.VoronoiIndex4[count*pnum:count*pnum + pnum]
    weights = num_coord(ql.num_connections, ydata, count, 1)
    plt.scatter(ql.num_connections, ydata, s=weights)
    plt.xlim(-1,13)
    plt.ylim(0,10)
    plt.xlabel('Q6-num-connections')
    plt.ylabel('VoronoiIndex4')
    plt.title("2Dplot,Q6,VI4-mc150-" + str(count))
    plt.savefig("Particle_dist-" + str(count) + ".png")    
    plt.clf()
    count += 1


#%%


for f in range(fnum):
    array = [0] * pnum
    for p in range(pnum):
        xdata = df.VoronoiIndex4[pnum*f+p]/df.Coordination[pnum*f+p]
        array[p] = xdata
    plt.hist(array, bins=20)
    plt.ylim(0,600)
    #plt.xlim(0,0.7)
    plt.xlabel('VI4-Coord-Ratio')
    plt.title("VI4-Coord-Ratio, frame " + str(count))
    plt.savefig("particle-dist-" + str(f) + ".png")    
    plt.clf()
    


#%%

#taking area as an array from a csv file, loading data

import csv
region = []
area = []
order = []

with open('VIarea.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for a in csv_reader:
        area.append([float(x) for x in a])
        

with open('VIregion.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for a in csv_reader:
        region.append([float(x) for x in a])


with open('VIorder.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for a in csv_reader:
        order.append([float(x) for x in a])
        

        
        

#%%

#Processing the arrays

partArray = [0] * pnum
totalArray = [0] * pnum

#make array with surface area of faces with four edges
#Additionally adding a threshhold for the surface area

for b in range(len(area)):
    if((int(order[b][0]) == 4) and area[b][0] <= 0.3):
        partArray[int(region[b][0])] += area[b][0]
        

for b in range(len(area)):
    totalArray[int(region[b][0])] += area[b][0]
    
fractionArray = [0] * pnum

for d in range(pnum):
    fractionArray[d] = partArray[d]/totalArray[d]

#%%
test =[]
for b in range(len(area)):
    if(int(order[b][0]) == 4): #and area[b][0] <= 0.4):
        test.append(area[b][0])
        
plt.plot(test)
plt.title("Area of 4-edged faces, frame 199")
plt.show()
        

#%%
#Plotting Surface area

mybins = np.linspace(0,0.3,60)
plt.hist(fractionArray,bins=mybins)
plt.xlabel("Surface area of VoronoiIndex4 to total Surface Area")
plt.ylabel("Num-particles")
plt.xlim(0, 0.3)
#plt.ylim(0,500)
plt.title("Frame 1750, relax484, 0.3 threshold")
plt.show()

#%%
#Graphing a 2d plot of the coordination and voronoiIndex4

frame = 1000

data = df.VoronoiIndex4[pnum*frame:pnum*(frame + 1) - 1]
coord = df.Coordination[pnum*frame:pnum*(frame + 1) - 1]
plt.plot(data, coord)
plt.show()

#%%

#Make string for particles with a threshold of surface area for VI4
string = ""
threshold = 0.001
for a in range(pnum):
    if(fractionArray[a] <= threshold):
        string += "ParticleIdentifier == " + str(a)
        if(a != pnum-1):
            string += " || "


#%%
#Identifying crystal liquid interface
import freud
import gsd.hoomd
import matplotlib.pyplot as plt
 

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
    ptype = [0] * pnum
    #counting the number of connections 
    for a in nlistdf.i:
        ptype[a] += 1
    #Removing partilces in bulk and liquid
    count = 0
    for b in range(len(nlistdf.j)):
        if(qlnc[nlistdf.j[b]] < 8): #Saying this neighbor particle is liquid
            pair[0] += 1
        else:
            pair[1] += 1 #Saying this neighbor particle is crystal
        #reset for each particle
        count += 1
        if(count == ptype[nlistdf.i[b]]):
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

#Making a gsd File based on particle positions


setnum = 150
particlenum = 2048
infile = '/Users/felixsong/Desktop/Ovito/mc55_2048_'+str(setnum)+'.gsd'
outfile = '/Users/felixsong/Desktop/Ovito/sim/cluster.gsd'
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
traj = gsd.hoomd.open(infile, 'rb')
framenum = len(traj)
trajout = gsd.hoomd.open(outfile, 'wb')
count = 0 
for framedata in (traj):
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    s = gsd.hoomd.Snapshot()
    s.particles.N = 2048 
    s.particles.position = framedata.particles.position
    s.configuration.box = [12.4927,12.4927,12.4927,0,0,0]
    s.particles.charge = interfaceP[count]
    count += 1
#     s.particles.mass = tempdf['cluster_idx78pp']
    trajout.append(s)

#%%

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


df1 = pd.DataFrame(datalist,columns=['size','rg','aveql','frame','particle','Particle-State','Position']) 


#%%



    

