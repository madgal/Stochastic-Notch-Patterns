import numpy as np
import matplotlib.pyplot as plt
from general import *
from datetime import datetime
import pandas as pd

odetype='ND'
lattice='2Dsquare'

## The parallel part start
from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

nAmpList=[]
nmax=100
step=2
if myrank==0 or myrank==2:
	nTypeL=['white']
	if nprocs>2:	
		if myrank==0:
			nAmpList = list(np.arange(0,110,10))
		elif myrank==2:
			nAmpList = list(np.arange(110,210,10))
	else:
		nAmpList = list(np.arange(0,210,10))
elif myrank==1 or myrank==3:
	nTypeL=['shot']
	if nprocs>2:	
		if myrank==1:
			nAmpList = list(np.arange(0,11,1))
		elif myrank==3:
			nAmpList = list(np.arange(11,21,1))
	else:
		nAmpList = list(np.arange(0,21,1))
else:
	exit()
#-----------------------------------
def add(arr,val):
	for el in arr:
		if el[0]==val[0] and el[1]==val[1]:
			return False
	return True

def flip(ics,rv,cv):
	rv = rv[0]
	cv = cv[0]
	maxICS=10000
	for key in ['I','N','D']:
		ics[key][0][rv][cv]=np.random.uniform(0,maxICS,1)[0]
		
	return ics
#-----------------------------------

############ Start simulation ###########
latsizeList=[16]
seedList=pd.read_csv(filepath_or_buffer='seedsForSims.txt',header=None).values[:,0]

simNum=1
tmax=5000
dt =0.1

misL = {0:[3,13,26],1:[64,128,192],2:[230,243,253]}

for size in latsizeList:
	nrow=size
	ncol=size
	nlayer=1
	for nType in ['shot']:
		for nAmp in [0]:
			for mistakes in misL[myrank]:#[3,13,26,64,128,192,230,243,253]:# 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99
				for seedInd in range(20):
					ss = seedList[seedInd]
					print size,nType,nAmp,mistakes,seedInd
                        
					[dat1,stat1]=run_simulation(odetype,nrow,ncol,nlayer,lattice,300,seed=ss,special='SP',
       							noiseType=nType,noiseAmp=nAmp,saveTraj=True,delta_t=dt,
       					                savefig=False,initialSys=False,isFlag=False,
       					                dimensionless=False,step4save=1,
							maxICSN=10000,maxICSD=10000,maxICSI=2000)
					ics = {'type':odetype}
					for k in dat1:
						ics[k] = dat1[k][-1]
                        
					np.random.seed(ss)
					v=[]
					while len(v)<mistakes:
						newEl = np.array([np.random.randint(0,nrow,1),np.random.randint(0,nrow,1)])
						t=add(v,newEl)
						if t:
							v+=[newEl]
					v = np.array(v)
					for k in range(len(v)):
						ics = flip(ics,v[k][0],v[k][1])

					[dat1,stat1]=run_simulation(odetype,nrow,ncol,nlayer,lattice,tmax,seed=ss,special=False,
       							noiseType=nType,noiseAmp=nAmp,saveTraj=True,delta_t=dt,
       					                savefig=False,initialSys=ics,isFlag=True,
       					                dimensionless=False,step4save=1,
							maxICSN=10000,maxICSD=10000,maxICSI=2000)
					for k in dat1:
						for i in range(len(dat1[k])):
							dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
					df_Final = pd.DataFrame(dat1)
					df_Final['run'] = seedInd
                        
					titleT = "data_ND_mcT_trajs_mistakes/traj_"+str(nrow)+"x"+str(ncol)+"x"+str(nlayer)+"_m"+str(mistakes)+"_s"+str(seedInd)+".dat"
					df_Final.to_csv(path_or_buf=titleT,index=False)
