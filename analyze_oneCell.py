## Written by Madeline Galbraith
## Last edited: July 2022

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib
import json

from aux_pseudopotential import *

tstart=100000
res=[]
Nstd,Dstd=[],[]
NextL,DextL=[],[]

count=0
for file in os.listdir("data_gil/oneCell/"):
    if 'json' in file:
        ## read in file as a JSON format
        filein = open('data_gil/oneCell/'+file)
        data = json.load(filein)
        
        ## select the D and N trajectories from gillespie one cell
        Dc = data['res']['D'][-tstart:]
        Nc = data['res']['N'][-tstart:]

        nxt = data['Next']
        dxt = data['Dext']
        prob=(np.mean(Dc)<np.mean(Nc))*1.
        nstd=np.std(Nc)
        dstd=np.std(Dc)
        if (nxt in NextL) and (dxt in DextL):
                indN = np.argwhere(nxt==np.array(NextL))[:,0]
                indD = np.argwhere(dxt==np.array(DextL))[:,0]
                ind = np.intersect1d(indD,indN)
                if len(ind)>0:
                    ind = ind[0]
                    res[ind]+=prob
                    Nstd[ind]+=[nstd]
                    Dstd[ind]+=[dstd]
                    #print 'A',prob,res[ind],Nstd[ind]
                else:
                    NextL+=[nxt]
                    DextL+=[dxt]
                    res+=[prob]
                    Nstd+=[[nstd]]
                    Dstd+=[[dstd]]
            
        else:
                NextL+=[nxt]
                DextL+=[dxt]                
                res+=[prob]
                Nstd+=[[nstd]]
                Dstd+=[[dstd]]
        count+=1
        if count%100==0:
            print count

#### Create a meshgrid so it can be plotted as a filled contour
x,y,z=[],[],[]
xun = np.unique(DextL)
for el in xun:
	inds = np.argwhere(el==DextL)[:,0]
	if len(inds)>0:
		tmpy,tmpz = np.array(NextL)[inds],np.array(res)[inds]
		total = np.array(Nstd)[inds]
		listx,listy,listz=[],[],[]
		for i in range(len(tmpy)):## no duplicates are ensured
			listx+=[el]
			listy+=[tmpy[i]]
			listz+=[tmpz[i]*1./len(total[i])]
		x+=[list(listx)]
		y+=[list(listy)]
		z+=[list(listz)]





print 'A'
print type(Dstd),type(Dstd[0]),np.array(Dstd).shape
Dstd = np.array(Dstd).flatten()
Nstd = np.array(Nstd).flatten()
print type(Dstd),type(Dstd[0]),np.array(Dstd).shape
stdHistD = np.histogram(Dstd,bins=50)
stdHistN = np.histogram(Nstd,bins=50)

print 'B'
results = {'Dstd':{'bins':list(stdHistD[1]),'freq':list(stdHistD[0])},'Nstd':{'bins':list(stdHistN[1]),'freq':list(stdHistN[0])},'Dext':x,'Next':y,'prob':z}

json_data =json.dumps(results,indent=4)
with open('gill1_analysis.json','w') as outfile:
	outfile.write(json_data)
