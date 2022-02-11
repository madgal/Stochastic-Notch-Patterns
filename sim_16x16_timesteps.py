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
if myrank==0:
	nTypeL=['white']
elif myrank==1:
	nTypeL=['shot']
else:
	print "Possible issue"
#	exit()
############ Start simulation ###########
latsizeList=[4,8,32]
seedList=pd.read_csv(filepath_or_buffer='seedsForSims.txt',header=None).values[:,0]

simNum=1
tmax=5000
dt =0.1
#for size in latsizeList:
for dt in [0.02]:#,0.01]:#,0.001]:
	size=16
	nrow=size
	ncol=size
	nlayer=1
	for nType in ['white']:#nTypeL:
		if nType=='shot':
			nAmpList =[0,7,13,20]
		elif nType=='white':
		    if dt==0.02:
			nAmpList =[26]
			#nAmpList =[14,26,40]# 700, 1300, 2000
		    elif dt==0.01:
			nAmpList =[7,13,20]# 700, 1300, 2000
		for nAmp in nAmpList:
			for seedInd in [0]:#range(1):
				ss = seedList[seedInd]
				[dat1,stat1]=run_simulation(odetype,nrow,ncol,nlayer,lattice,tmax,seed=ss,special=False,
       							noiseType=nType,noiseAmp=nAmp,saveTraj=True,delta_t=dt,
       					                savefig=False,initialSys=False,isFlag=False,
       					                dimensionless=False,step4save=1,
							maxICSN=10000,maxICSD=10000,maxICSI=2000)
				for k in dat1:
					for i in range(len(dat1[k])):
						dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
				df_Final = pd.DataFrame(dat1)
				df_Final['run'] = seedInd

				if dt==0.02:
					titleT = "data_ND_mcT_trajs_time/traj_"+str(dt)+"_"+str(nrow)+"x"+str(ncol)+"x"+str(nlayer)+"_"+nType+"_n"+str(nAmp*5)+".dat"
				elif dt==0.01:
					titleT = "data_ND_mcT_trajs_time/traj_"+str(dt)+"_"+str(nrow)+"x"+str(ncol)+"x"+str(nlayer)+"_"+nType+"_n"+str(nAmp*10)+".dat"
				df_Final.to_csv(path_or_buf=titleT,index=False)
