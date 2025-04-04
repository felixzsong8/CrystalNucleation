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
#This is for PTM
#clist = ['x','y','z','StructureType', 'id']
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
    if(df.VoronoiIndex4[a] >= 2) and (df.Coordination[a]>=14):
        index4[df.frame[a]] += 1
        


plt.plot(index4)
plt.title('Nucleation Sim: VoronoiIndex4 >= 2 and Coordination >= 14')
plt.xlabel("Frame")
plt.ylabel("Number of Crystal Particles")
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

infile = "/Users/felixsong/Desktop/Ovito/mc55_2048_150.gsd" # your file path 
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
#Function to create number of points at a certain coordinate
#Where the ydata is using df

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

#Same function but with no df variable
#C is the multiplier for the wieghts
def num_coord_nodf(xdata,ydata,count,c):
    dict = {}
    array = [0] * len(xdata)

    for i in range(len(xdata)):
        if (xdata[i],ydata[i]) in dict:
            dict[xdata[i],ydata[i]] += 1
        else:
            dict[xdata[i],ydata[i]] = 1

    for j in range(len(xdata)):
        array[j] = c*dict[xdata[j],ydata[j]]
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
framenum = 199

with open('/Users/felixsong/Desktop/Ovito/areaData/mc150/VIarea' + str(framenum) +'.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for a in csv_reader:
        area.append([float(x) for x in a])
        

with open('/Users/felixsong/Desktop/Ovito/areaData/mc150/VIregion' + str(framenum) +'.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for a in csv_reader:
        region.append([float(x) for x in a])


with open('/Users/felixsong/Desktop/Ovito/areaData/mc150/VIorder' + str(framenum) +'.csv') as csv_file:
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
    if((int(order[b][0]) == 4) and area[b][0] <= 1):
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
        
plt.hist(test, bins = 20)
plt.title("Area of 4-edged faces, frame 199")
plt.show()
        

#%%
#Plotting Surface area
fn = 199
mybins = np.linspace(0,0.05,60)
plt.hist(fractionArray[fn*pnum: fn*pnum +pnum-1],bins=mybins)
plt.xlabel("Surface area of VoronoiIndex4 to total Surface Area")
plt.ylabel("Num-particles")
plt.xlim(0, 0.05)
plt.ylim(0,600)
plt.title("Frame 199, mc150, 0.2 threshold")
plt.show()

#%%
#Graphing a 2d plot of the coordination and voronoiIndex4

frame = 1000

data = df.VoronoiIndex4[pnum*frame:pnum*(frame + 1) - 1]
coord = df.Coordination[pnum*frame:pnum*(frame + 1) - 1]
plt.plot(data, coord)
plt.show()

#%%
#Make a array for the fraction of each cells VI4
import csv

fractionArray = []
for a in range(fnum):
    with open('/Users/felixsong/Desktop/Ovito/areaData/mc150-0.12/FractionArray' + str(a) + '.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        frac = []
        for b in csv_reader:
            frac.append([float(x) for x in b])
        for c in frac:
            fractionArray.append(c[0])

#%%

#Make string for particles with a threshold of surface area for VI4
string = ""
facelist = []
threshold = 0.05
for a in range(pnum):
    if(ydata[a] <= threshold and ydata[a] != 0):
        string += "ParticleIdentifier == " + str(a)
        facelist.append(a)
        if(a != pnum-1):
            string += " || "


#%%
#Identifying crystal liquid interface
#
#
#
#
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
#Recursive function that, given a starting particle, will traverse through all neighbor cystal particles
#and put them into an array, stops at liquid interface

def clusterTrav(particleStart, visitedParticles, interfaceP, nlist, frame):
    if particleStart in visitedParticles or interfaceP[frame][particleStart] == 3:
        return []
    if interfaceP[frame][particleStart] == 1: #Once it reaches this particle
        visitedParticles.append(particleStart)
        return []
    visitedParticles.append(particleStart)
    for a in nlist[particleStart]:
        clusterTrav(a, visitedParticles, interfaceP, nlist, frame)
    return visitedParticles
        
        
    


#%%
#update: adding cluster numbering and another property to find which cluster they are next to

ClusterArray = []
frame = 0
for framedata in traj:
    index = 0 #to id each cluster
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    qlnc = ql.num_connections
    nnlist = []
    for i,j in ql.nlist:
        nnlist.append([i,j])
    nlistdf = pd.DataFrame(nnlist,columns=['i','j'])
    clusterIndex = [-1] * pnum
    ptype = [0] * pnum
    nlist = []
    array = []
    count = 0
    #counting the number of connections and putting it into array
    for a in nlistdf.i:
        ptype[a] += 1
    for b in nlistdf.j:
        if(count == ptype[nlistdf.i[b]]):
            nlist.append(array)
            count = 0
            array = []
        array.append(b)
        count += 1
    nlist.append(array)
    #go through each particle and run recursive function if they aren't labeled
    for c in range(pnum):
        if(clusterIndex[c] == -1):
            parray = clusterTrav(c, [], interfaceP, nlist, frame)
            if(len(parray) != 0):
                for d in parray:
                    if(clusterIndex[d] != -1 and interfaceP[frame][d] == 2):
                        clusterIndex[d] = clusterIndex[d]*1000 + index
                    else:
                        clusterIndex[d] = index
                index += 1
    ClusterArray.append(clusterIndex)
    frame += 1
    print(frame)
    
#%%
string = ""
array = p
for a in array:
    string += "ParticleIndex == " + str(a)
    if(a != pnum-1):
        string += " || "


    
#%%
#Making array of cluster Id based on surface interface

def check(c, nlist, array):
    for a in nlist:
        if(qlnc[a] >= 8 and qlci[a] != qlci[c]):
            return array[a]
    return qlci[c]
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
#Here, we can see in frame 69 that one particle is not clustered correctly in qlci
#We can see that 643 and 131 are both solid and nearest neighbors but not in same cluster
count = 0
framedata = traj[34]
ql.compute(system=framedata, neighbors={'num_neighbors': 12})
qlnc = ql.num_connections
qlci = ql.cluster_idx
nnlist = []
for i,j in ql.nlist:
    nnlist.append([i,j])
nlistdf = pd.DataFrame(nnlist,columns=['i','j'])
for a in range(len(qlci)):
    if(clusterId[34][a] == 659):
        count = a
#%%
for c in range(12):
    print(nlistdf.j[12*211+c])
   # print(framedata.particles.charge[nlistdf.j[12*473+c]])
    print(qlnc[nlistdf.j[12*211+c]])

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
    s.particles.charge = clusterId[count*pnum:(count+1)*pnum]
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

#PLotting charge and Surface area
import csv
infile = "/Users/felixsong/Desktop/Ovito/sim/cluster.gsd" # your file path 
traj = gsd.hoomd.open(infile, 'rb')

fig = plt.figure()
framenum = len(traj)
count = 0
for framedata in traj:
    xdata = framedata.particles.charge
    frac = []
    with open('/Users/felixsong/Desktop/Ovito/areaData/mc150/FractionArray' + str(count) + '.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for a in csv_reader:
            frac.append([float(x) for x in a])
    ydata = []
    for a in frac:
        ydata.append(a[0])
    weights = num_coord_nodf(xdata, ydata, count, 5)
    plt.scatter(xdata, ydata ,s=weights)
    plt.xlim(-1,4)
    plt.ylim(0,0.3)
    plt.xlabel('Particle interface')
    plt.ylabel('SA Fraction of VI4 w/ 0.3 filter')
    plt.title("Particles at the interface vs VI4")
    plt.savefig("Particle_dist-" + str(count) + ".png")    
    plt.clf()
    count += 1


#%%
#Using our threshold parameters to analyze particles


threshold = 0.01
framenum = 0
tw = [0]*fnum 
pl = []
for f in range(fnum):
    partlist = []
    data = []
    xdata = []
    ydata = []
    for particles in range(pnum):
        if(fractionArray[framenum*pnum + particles] <= threshold and fractionArray[framenum*pnum + particles] != 0):
            partlist.append(particles)
    for a in partlist:
        xdata.append(df.Coordination[framenum*pnum + a])
        ydata.append(df.VoronoiIndex4[framenum*pnum + a])
        if(df.Coordination[framenum*pnum + a] >= 12 and df.VoronoiIndex4[framenum*pnum + a] >= 2):
            tw[f] += 1
            data.append(a)
    weights = num_coord_nodf(xdata, ydata, framenum, 3)
    plt.scatter(xdata, ydata ,s=weights)
    plt.xlim(10,20)
    plt.ylim(-1,7)
    plt.xlabel('Coordination')
    plt.ylabel('VoronoiIndex4')
    plt.title("VI4 and Coordination of cells with " + str(threshold) + " on faces and 0.2 threshold on fraction Index")
    plt.savefig("Particle_dist-" + str(framenum) + ".png")    
    plt.clf()
    pl.append(data)
        
    framenum += 1
    
#%%

#Plotting the number of cells with parameters defined above
#Threshold of faces and fraction
 
plt.plot(tw)
plt.title("Particles identified as crystals w/threshold at 0.08-0.2")
plt.show()


#%%
#return a string that has all of the particles with certain parameters

string = ""
stringnot = ""
facelist = []
threshold = 0.01
framenum = 0
for a in range(pnum):
    if(fractionArray[framenum*pnum + a] <= threshold and fractionArray[framenum*pnum + a] != 0.0 and df.Coordination[framenum*pnum + a] >= 12 and df.VoronoiIndex4[framenum*pnum + a] >= 2):
        string += "ParticleIdentifier == " + str(a + 1)
        facelist.append(a)
        if(a != pnum-1):
            string += " || "
    else:
        stringnot += "ParticleIdentifier == " + str(a + 1)
        if(a != pnum-1):
            stringnot += " || "


#%%

#Make a array for the fraction of each cells VI4
import csv
point = 0.05
for a in range(11):
    point = round(point,2)
    fractionArray = []
    for a in range(fnum):
        with open('/Users/felixsong/Desktop/Ovito/areaData/mc150-' + str(round(point,2)) + '/FractionArray' + str(a) + '.csv') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            frac = []
            for b in csv_reader:
                frac.append([float(x) for x in b])
            for c in frac:
                fractionArray.append(c[0])
    threshold = 0.005
    for b in range(10):
        framenum = 0
        tw = [0]*fnum 
        for f in range(fnum):
            partlist = []
            xdata = []
            ydata = []
            for particles in range(pnum):
                if(fractionArray[framenum*pnum + particles] <= threshold and fractionArray[framenum*pnum + particles] != 0):
                    partlist.append(particles)
            for a in partlist:
                xdata.append(df.Coordination[framenum*pnum + a])
                ydata.append(df.VoronoiIndex4[framenum*pnum + a])
                if(df.Coordination[framenum*pnum + a] >= 12 and df.VoronoiIndex4[framenum*pnum + a] >= 2):
                    tw[f] += 1
            framenum += 1
        plt.plot(tw)
        #plt.xlim(10,20)
        #plt.ylim(-1,7)
        plt.xlabel('Frames')
        plt.ylabel('Particles Identified')
        plt.title("Crystal Identification with face-" + str(round(point,2)) + "-frac-" + str(round(threshold,3)))
        plt.savefig("/Users/felixsong/Desktop/Ovito/plots/mc150-n0/CrystalIdentification-" + str(round(point,2)) + "-" + str(round(threshold,3)) + ".png")    
        plt.clf()
        threshold += 0.005
    point += 0.01
    
    
#%%
#First function reads a line and splits it into an array of nums
#Second function just defines whether number is between a range
def readNum(line):
    text = []
    string = ""
    for length in range(len(line)):
        if(line[length] != " "):
            string += line[length]
        else:
            text.append(float(string))
            string = ""
    return text

def inRange(num):
    if(num >= 85 and num <= 100):
        return True
   # elif(num >= 125 and num <= 145):
    #    return True
    elif(num >= 160 and num <= 180):
        return True
    else:
        return False


#%%
#Takes a text file of Bond angles, Tiplets, bond pairs, and Voronoi Order
#and forms into array, particleTrack, with shows every particle in every frame
#and the number of angles within upper and lower that go through 4 sided faces

numf = 200
partnum = 2048

particleTrack = [0] * partnum*numf
for a in range(numf):
    count = 0
    with open('test.' + str(a) + '.txt') as f:
        for header in range(3):
            lines = f.readline()
        numList = readNum(lines)
        for b in range(partnum):
            switch = True
            while(switch):
                if(len(numList) == 0):
                    switch = False
                elif(numList[5] != count):
                    count += 1
                    switch = False
                else:
                    if(numList[1] + numList[2] == 0 and inRange(numList[3])): 
                        particleTrack[partnum*a + b] += 1 
                    lines = f.readline()
                    numList = readNum(lines)
    
                    
    
#%%
#Puts into array the number per frame of particles with connections

track = [0] * numf
trackLow = [0] * numf

for g in range(numf):
    for h in range(partnum):
        if(particleTrack[partnum*g + h] >= 1):
            track[g] += 1
        if(particleTrack[partnum*g + h] >= 2):
            trackLow[g] += 1

    
#%%

plt.plot(track)
plt.xlabel("frames/time") 
plt.ylabel("num particles")
plt.title("Number of particles with bond angles in range per frame \n range: 80-105, 125-145, 160-180, threshold-0.3 \n At least one bond angle")
plt.show()

plt.plot(trackLow)
plt.xlabel("frames/time") 
plt.ylabel("num particles")
plt.title("Number of particles with bond angles in range per frame \n range: 80-105, 125-145, 160-180, threshold-0.3 \n At least two bond angle")
plt.show()

#%%

#Convert particleTrack to array of 200 * 2048
data = []
for i in range(numf):
    slot = []
    for j in range(partnum):
        if(particleTrack[partnum*i + j] >= 2):
            slot.append(1) 
        else: 
            slot.append(0)
    data.append(slot)

#%%

#Make a new GSD file with an included property for BondAngle
import freud
import gsd.hoomd
import matplotlib.pyplot as plt

particlenum = 2048
infile = '/Users/felixsong/Desktop/Ovito/mc55_2048_150.gsd'
outfile = '/Users/felixsong/Desktop/Ovito/sim/BondAng2.gsd' 
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
    s.particles.charge = data[count]
    s.particles.Q6 = interfaceP[count]
    count += 1
    trajout.append(s)
    
#%%
    
#Comparing PTM with our parameters
#Read data file where structuretype is ptm


filename = "/Users/felixsong/Desktop/Ovito/mc150_ptm"
#This is for PTM
clist = ['x','y','z','StructureType', 'id']
df, pnum, fnum = readOvito(filename, clist)
arrayPTM = [0] * fnum
arrayRel = [0] * fnum
for i in range(fnum):
    for j in range(pnum):
        if(df.StructureType[pnum*i + j] == 1):
            arrayPTM[i] += 1
    arrayRel[i] = arrayPTM[i]/track[i]

#%%

#load Q6
import freud
import gsd.hoomd
import matplotlib.pyplot as plt
infile = "/Users/felixsong/Desktop/Ovito/sim/cluster.gsd" # your file path 
traj = gsd.hoomd.open(infile, 'rb')
arrayQ6 = [0] * fnum
arrayQRel = [0] * fnum
ql = freud.order.SolidLiquid(l=6,q_threshold=0.7,solid_threshold=8)
frame = 0
for framedata in traj:
    count = 0
    ql.compute(system=framedata, neighbors={'num_neighbors': 12})
    qlnc = ql.num_connections
    for i in qlnc:
        if(i >= 8):
            count += 1
    arrayQ6[frame] += count
    arrayQRel[frame] = arrayQ6[frame]/track[frame]
    frame += 1
    



#%%


#Plotting 3D plot between PTM, VI, frame
from mpl_toolkits import mplot3d
fig = plt.figure()

ax = plt.axes(projection = '3d')
xdata = arrayQ6
ydata = track
zdata = range(200)
ax.scatter3D(xdata, ydata, zdata, c = zdata, s=3)
plt.xlim(0,2048)
plt.ylim(0,2048)
ax.set_zlim(0,200)
ax.set(xlabel='Q6', ylabel='Bond-Angle-Parameter',zlabel='Frame')
ax.set_title("3Dplot,mc150")
plt.savefig("2plot.png")    
plt.clf()

#%%

#2D graph showing the correlation between PTM and BondAngle
#ArrayQ6 is Q6/VI

xdata = range(200)
ydata1 = track
ydata2 = arrayQ6
ydata3 = arrayPTM
ydata4 = trackLow
plt.plot(xdata, ydata1, label='BondAng-High')
plt.plot(xdata, ydata2, label='Q6')
plt.plot(xdata, ydata3, label='PTM')
plt.plot(xdata, ydata4, label='BondAng-Low')
leg = plt.legend(loc='upper left')
#plt.xlim(0,2048)
#plt.ylim(0,2048)

plt.xlabel('Frame')
plt.ylabel('Ratio of Q6 to BondAng')
plt.title("mc150, ") 
plt.savefig("2plot.png")    
plt.clf()
    
    
    
    
    
    