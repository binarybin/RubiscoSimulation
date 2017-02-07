from __future__ import division
from sys import argv
import numpy as np
from numpy import log
np.set_printoptions(threshold=np.nan)
from glob import glob



# In[4]:

def findBegin(lines):
    for idx, line in enumerate(lines):
        if line == '[Polymerlist]\n':
            return idx

def findEnd(lines):
    for idx, line in enumerate(lines): 
        if line == 'END\n':
            return idx

def splitSimSumo(lines):
    for idx, line in enumerate(lines):
        if line[:4] == "Sumo":
            return lines[:idx], lines[idx:]
            
def parseHead(lines, size):
    layer = np.zeros((size,size), int)
    for idx, line in enumerate(lines):
        if line[0] != "S":
            headString = line.split(" ")[0]
            x = int(headString.split(",")[0][1:])
            y = int(headString.split(",")[1][:-1])
            layer[x][y] = 1
    return layer

def get_polys(lines):
    for idx, line in enumerate(lines):
        if line == 'Two dimensional single layer space\n':
            start = idx+1
        if line == 'No bond info for this geometry\n':
            end = idx - 3
            break

    polylines = lines[start:end]
    
    polyodd  = [pline for idx, pline in enumerate(polylines) if idx %2 == 0]
    polyeven = [pline for idx, pline in enumerate(polylines) if idx %2 != 0]
    
    polyinfo = zip(polyodd, polyeven)
    
    sims = {}
    sumos = {}
    simcount = 1
    sumocount = 1
    for element in polyinfo:
        identifier = element[0].split()[0]
        points = element[1].split()
        pointinfo = [info.split('(')[1].split(')')[0].split(',') for info in points]
        pointfull = [(int(info[0]), int(info[1])) for info in pointinfo]
        if identifier == 'Sim':
            sims[simcount]=pointfull
            simcount += 1
        elif identifier == 'Sumo':
            sumos[sumocount]  = pointfull
            sumocount += 1
        else:
            raise Exception("Something wrong with the type")
    return sims, sumos

def get_layers(lines, sims, sumos):
    for idx, line in enumerate(lines):
        if line == 'END\n':
            size = idx - 3
            break
    simlayer = np.zeros((size,size), int)
    sumolayer = np.zeros((size,size), int)
    for simkey in sims.keys():
        for pt in sims[simkey]:
            x, y = pt
            simlayer[x,y] = simkey
    for sumokey in sumos.keys():
        for pt in sumos[sumokey]:
            x, y = pt
            sumolayer[x,y] = sumokey

    return simlayer, sumolayer


def classify_clusters(sims, sumos, siml, sumol):
    simclusters = []
    sumoclusters = []
    simset = set(sims.keys())
    sumoset = set(sumos.keys())
    while sumoset:
        sumo = sumoset.pop()
        tasklist = [(sumo, 'o')]
        localsims = set()
        localsumos = set([sumo])
        for task in tasklist:
            if task[1] == 'o': # that is a sumo
                pts = sumos[task[0]]
                for pt in pts:
                    x, y = pt
                    simid = siml[x,y]
                    if simid not in localsims and simid != 0:
                        tasklist.append((simid, 'i'))
                        localsims.add(simid)
            else: # that is a sim
                pts = sims[task[0]]
                for pt in pts:
                    x, y = pt
                    sumoid = sumol[x,y]
                    if sumoid not in localsumos and sumoid != 0:
                        tasklist.append((sumoid, 'o'))
                        localsumos.add(sumoid)
        simset = simset - set(localsims)
        sumoset = sumoset - set(localsumos)
        sumoclusters.append(localsumos)
        simclusters.append(localsims)
    #    print sumoclusters
    return sumoclusters

def analyze(lines):
    sims, sumos = get_polys(lines)
    siml, sumol = get_layers(lines, sims, sumos)
    sumoclusters = classify_clusters(sims, sumos, siml, sumol)
    sizes = np.array([len(ele) for ele in sumoclusters])
    return sizes



size = 50
f = open(argv[1], 'r')
line = f.readlines()
propraw = argv[1]
nepyc = propraw.split("nsim")[1].split('_')[0]
nrubi = propraw.split("nsumo")[1].split('_')[0]
lepyc = propraw.split("lsim")[1].split('_')[0]
lrubi = propraw.split("lsumo")[1].split('_')[0]
beta = propraw.split("beta_")[1].split('_')[0]
gamma = propraw.split("gamma_")[1].split(".txt")[0]
step = propraw.split("step_")[1].split("_")[0]
sizes = analyze(line)

#print sizes
fwrite = open('sizeDistri.txt', 'a')
sizestr = "["
for num in sizes:
    sizestr += str(num)
    sizestr += " "
sizestr += "]"
fwrite.write(nepyc+"\t"+nrubi+"\t"+lepyc+"\t"+lrubi+"\t"+beta+"\t"+gamma+"\t"+step+"\t SIZES: "+sizestr+"\n")
fwrite.close()
