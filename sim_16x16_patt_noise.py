import numpy as np
import matplotlib.pyplot as plt
from general import *
from datetime import datetime
import time
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
			#nAmpList = [0,50]#list(np.arange(0,110,10))
			nAmpList=[10,20,30,40,60,70,80,90]
		elif myrank==2:
			#nAmpList = [100,150,200]#list(np.arange(110,210,10))
			nAmpList=[110,120,130,140,160,170,180,190]
	else:
		nAmpList = [0,50,100,150,200]#list(np.arange(0,210,10))
elif myrank==1 or myrank==3:
	nTypeL=['shot']
	if nprocs>2:	
		if myrank==1:
			#nAmpList = [0,5]#list(np.arange(0,11,1))
			nAmpList=[1,2,3,4,6,7,8,9]
		elif myrank==3:
			#nAmpList = [10,15,20]#list(np.arange(11,21,1))
			nAmpList=[11,12,13,14,16,17,18,19]
	else:
		nAmpList = [0,5,10,15,20]#list(np.arange(0,21,1))
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

recv ={'N':np.max(dat0['N'][-1]),'D':np.min(dat0['D'][-1]),'I':np.max(dat0['I'][-1])}
senv ={'N':np.min(dat0['N'][-1]),'D':np.max(dat0['D'][-1]),'I':np.min(dat0['I'][-1])}

shapeV = dat0['N'][-1].shape

print recv,senv,'\n'

def getLatticeICS(pattern):
	
	ics={}
	ics['N']=np.zeros(shapeV)
	ics['I']=np.zeros(shapeV)
	ics['D']=np.zeros(shapeV)

	if pattern==0:
		## only one corner is receivers
		for k in ics:
			ics[k][0,:8,:8]=recv[k]	
			ics[k][0,:8,8:]=senv[k]	
			ics[k][0,8:]=senv[k]	
	elif pattern==1:
		## only one corner is senders
		for k in ics:
			ics[k][0,:8,:8]=senv[k]	
			ics[k][0,:8,8:]=recv[k]	
			ics[k][0,8:]=recv[k]	
	elif pattern==2:
		## all senders
		for k in ics:
			ics[k][0,:,:]=senv[k]	
	elif pattern==3:
		## all receivers
		for k in ics:
			ics[k][0,:,:]=recv[k]	
	elif pattern==4:
		## one line receivers
		for k in ics:
			ics[k][0,:,:]=senv[k]	
			ics[k][0,2,:]=recv[k]	
	elif pattern==5:
		## one line senders
		for k in ics:
			ics[k][0,:,:]=recv[k]	
			ics[k][0,2,:]=senv[k]	
	elif pattern==6:
		## half senders
		for k in ics:
			ics[k][0,:8,:]=senv[k]	
			ics[k][0,8:,:]=recv[k]	
	elif pattern==7:
		## nucleating site
		randx,randy = np.random.randint(0,16,2)
		for k in ics:
			ics[k][0,:,:]=recv[k]	
			ics[k][0,randx,randy]=senv[k]	

	elif pattern==8:
		for k in ics:	
			ics[k]=dat0[k][-1]
				
	ics['type']='ND'

	return ics

	
ptitle={0:"cornerR",1:"cornerS",2:"S",3:"R",4:"lineR",5:"lineS",6:"Half"}

for size in latsizeList:
	nrow=size
	ncol=size
	nlayer=1
	for nType in nTypeL:
		for nAmp in nAmpList:
			for pattern in [8]:#[0,4,7,8]:#range(8):
				print size,nType,nAmp,pattern
				if pattern==7 or pattern==8:
					listL =  20
				else:
					listL=1#10
				for seedInd in range(listL):
					t0=time.time()
					ss = seedList[seedInd]
					ics=getLatticeICS(pattern)
					#print ics['N'],'\n'
                       
					[dat1,stat1]=run_simulation(odetype,nrow,ncol,nlayer,lattice,tmax,seed=ss,special=False,
       					noiseType=nType,noiseAmp=nAmp,saveTraj=True,delta_t=dt,
       				        savefig=False,initialSys=ics,isFlag=True,
       				        dimensionless=False,step4save=1,
					maxICSN=10000,maxICSD=10000,maxICSI=2000)
					for k in dat1:
						for i in range(len(dat1[k])):
							dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
					#if seedInd==0:
					df_Final = pd.DataFrame(dat1)
					df_Final['run'] = seedInd

					print time.time()-t0
					titleT = "data_ND_mcT_trajs_mar2021/patt/traj_"+str(nrow)+"x"+str(ncol)+"x"+str(nlayer)+"_"+str(nType)+"_n"+str(nAmp)+"_p"+str(pattern)+"_s"+str(seedInd)+".dat"
					df_Final.to_csv(path_or_buf=titleT,index=False)#,mode='a')
