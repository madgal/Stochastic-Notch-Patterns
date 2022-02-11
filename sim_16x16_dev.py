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
	for key in ics:
		ics[key][0][rv][cv]=np.random.uniform(0,maxICS,1)[0]
	return ics
#-----------------------------------

############ Start simulation ###########
latsizeList=[16]
seedList=pd.read_csv(filepath_or_buffer='seedsForSims.txt',header=None).values[:,0]

simNum=1
tmax=5000
dt =0.1
nrow,ncol=16,16
nlayer=1
[dat0,stat1]=run_simulation(odetype,nrow,ncol,nlayer,lattice,300,seed=13141,special='SP',
		noiseType='white',noiseAmp=0,saveTraj=True,delta_t=dt,
                savefig=False,initialSys=False,isFlag=False,
                dimensionless=False,step4save=1,
		maxICSN=10000,maxICSD=10000,maxICSI=2000)
dList={0:[0,1,2,3,4],1:[5,10,15,25],2:[50,75,100,200,500],3:[1000,1500,2000,5000]}
for size in latsizeList:
	nrow=size
	ncol=size
	nlayer=1
	for nType in ['shot']:
		for nAmp in [0]:
			for sdev in dList[myrank]:#[0,1,2,3,4,5,10,15,25,50,75,100,200,500,1000,1500,2000,5000]:
				for seedInd in range(20):
					ss = seedList[seedInd]
					print size,nType,nAmp,sdev,seedInd
                        
					np.random.seed(ss)
					ics = {'type':odetype}
					for k in dat0:
						ics[k] = dat0[k][-1]
						ics[k]=np.random.normal(ics[k],sdev,ics[k].shape)
						ics[k] = (ics[k]>=0)*ics[k]
                        
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
                        
					titleT = "data_ND_mcT_trajs_dev/traj_"+str(nrow)+"x"+str(ncol)+"x"+str(nlayer)+"_d"+str(sdev)+"_s"+str(seedInd)+".dat"
					df_Final.to_csv(path_or_buf=titleT,index=False)
