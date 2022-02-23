import numpy as np
import pandas as pd

## Allow the calculation to be parallelized
from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()


###################
###################
###################
## Parameters for notch system
N0 = 5.0e+2     ##molecules/hour 
D0 = 1.0e+3     ##molecules/hour 
I0 = 2.0e+2     ##molecules

kc = 5.0e-4     #1/hour/molecule
kt = 5.0e-5     #1/hour/molecule

nI = 2.0e+0        #hill coeff 

lambdaN = 2.0e+0   # no dimension
lambdaD = 0.0e+0   # no dimension

gamma = 1.0e-1   #1/hour
gamma_i = 5.0e-1 #1/hour

###################
###################
###################
def generate_lattice(nrow,ncol,nlayer,special,maxICS):
    ### Create the initial conditions for 
    ###         -'1D'
    ###    2D square lattice

    ### nlayer >1 means a 2d system
    ### nrow and ncol is the length and width of system
    nlayer=1
    a= np.zeros((nlayer,nrow,ncol))

    if not special:
        a= np.random.uniform(0,maxICS,[nlayer,nrow,ncol])
    elif special=='SP':
        a[nlayer/2][nrow/2][ncol/2]=3000
    
    return a
###################
###################
###################
def repackage(total):
	## rebuild dictionary so its easier to parse for 
	
	test = total[0]
	keys = test.keys()

	update={}
	for k in keys:
		update[k]=[]

	for i in range(len(total)):
		for k in keys:
			update[k]+=[total[i][k]]		

	return update
###################
###################
###################

def shiftedHill(X,X0,nX,lambdY):
	return lambdY+(1.0 - lambdY)/(1.0+(X/X0)**nX)
###################
###################
###################
def ts_update(dt,X,tempX,noiseType,noiseX):

	if noiseType=='shot' and noiseX!=0:
		sigma=np.sqrt(dt*X)*np.random.normal(0,noiseX,X.shape)
	elif noiseType=='white' and noiseX!=0:
		sigma=np.random.normal(0,noiseX,X.shape)*dt
	else:
		sigma=0.

	X = X+dt*tempX+sigma

	decX = np.sign(X)+1 ## decX is now 0 if X neg, 1 if X=0, and 2 if X pos
	X = X*decX/2. ## zeros out neg X vals, doesn't modif if X==0 or X is positive

	return X
###################
###################
###################
def getExt(X,lattice):
	if '1D' in lattice:
		## the value of Next and Dext is equal to the value of N or D in the neighbor 
		if X.shape[1]>1:
			## if the array has two columns
			return np.roll(X,1,axis=1)
		elif X.shape[2]>1:
			## if the array has two rows
			return np.roll(X,1,axis=2)
	elif lattice=='2Dsquare':
		## the value of Next and Dext is equal to the 1/4 of N or D in the neighbor as each cell has 4 neighbors
	    return 1./4.*(np.roll(X,1,axis=1)+np.roll(X,-1,axis=1)+np.roll(X,1,axis=2)+
			  np.roll(X,-1,axis=2))
###################
###################
###################
def run_simulation(tmax,delta_t,lattice,latticesize,noiseType,noiseAmp,special=None,ICS=None,step4save=None):
	if not step4save:
		step4save=1
	if not ICS:
		## get the initial conditions
		## if special is not set than conditions are randomized
		Ni=generate_lattice(latticesize,latticesize,1,special,10000)
		Di=generate_lattice(latticesize,latticesize,1,special,10000)
		Ii=generate_lattice(latticesize,latticesize,1,special,2000)

		varbls = {'N':Ni,'D':Di,'I':Ii}
	else:
		varbls=ICS

	total=[]
	for i in range(int(tmax/delta_t)):
		## Nc,Dc, and Ic are the Notch, Delta, and NICD values at time=t
		Nc,Dc,Ic = varbls['N'],varbls['D'],varbls['I']

		## calculate the value of Next and Dext, depends on whether the system is 2 cell or multicell
                Next = getExt(Nc,lattice)
                Dext = getExt(Dc,lattice)

		## temp_X =dX/dt, must multiply by dt to get the dX which will be added to value at time=t 
		temp_N = N0*shiftedHill(Ic,I0,nI,lambdaN) -Nc*(kc*Dc+gamma+kt*Dext) 
		temp_D = D0*shiftedHill(Ic,I0,nI,lambdaD) -Dc*(kc*Nc+gamma+kt*Next)
		temp_I = kt*Nc*Dext-gamma_i*Ic     

		## get value of N,D, and I at time=t+dt
		## also include the noise calculation
		Nc = ts_update(delta_t,Nc,temp_N,noiseType,noiseAmp)
		Dc = ts_update(delta_t,Dc,temp_D,noiseType,noiseAmp)
		Ic = ts_update(delta_t,Ic,temp_I,noiseType,0)## use noiseI=0 because we don't add noise here

		varbls={'N':Nc,'D':Dc,'I':Ic}
		if i%step4save==0:
			total.append(varbls)

	total = repackage(total)

	return total
###################
###################
###################
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
def getLatticeICS(pattern,recv,senv,shapeV):
	
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
				
	return ics
#-----------------------------------



	
##########################################
############ Start simulation ###########
##########################################
## Code will work with or without mpi
## setup to run the different types of noise on different processes
nAmpList={'white':{},'shot':{}}
nTypeL={}
if nprocs==1:
	nTypeL[0]=['white','shot']
	nAmpList['white'][0] = list(np.arange(0,2100,100))
	nAmpList['shot'][0] = list(np.arange(0,21,1))
elif nprocs<=3:
	nTypeL[0]=['white']
	nTypeL[1]=['shot']
	nAmpList['white'][0] = list(np.arange(0,2100,100))
	nAmpList['shot'][1] = list(np.arange(0,21,1))
elif nprocs>=4:
	nTypeL[0]=['white']
	nTypeL[1]=['white']
	nTypeL[2]=['shot']
	nTypeL[3]=['shot']
	nAmpList['white'][0] = list(np.arange(0,1100,100))
	nAmpList['white'][1] = list(np.arange(1100,2100,100))
	nAmpList['shot'][2] = list(np.arange(0,11,1))
	nAmpList['shot'][3] = list(np.arange(11,21,1))

if myrank not in nTypeL.keys():
	## if the process won't do anything than exit to save computation time
	exit()

## get seeds so data is consistent across runs
seedList=pd.read_csv(filepath_or_buffer='seedsForSims.txt',header=None).values[:,0]

## parameters for simulation
mistakeList=[13,64,230]
sdevList=[50,1500,5000]
pattList=[0,7]#0-7 for patterns listed above]
latsizeList=[16]
tmax=10000
dt =0.1
numberOfSims=20
lattice='2Dsquare'


## Set these flags to determine which simulation(s) you want to run
genMistakes=False
genDeviations=False
genPatterns=False
genRandom=True

for latticesize in latsizeList:
    for nType in nTypeL[myrank]:
	for nAmp in nAmpList[nType][myrank]:

	    if genMistakes:
	    	for mistakes in mistakeList:   
			for seedInd in range(numberOfSims):
				ss = seedList[seedInd]
				np.random.seed(ss)
				nrow=latticesize
				ncol=latticesize

				## Initial condition: Salt and pepper with discrete perturbations (mistakes in the lattice)
				## First generate the salt and pepper pattern
				datsp=run_simulation(tmax,dt,lattice,latticesize,nType,0,special='SP')

				##
				ics={}
				for k in datsp:
					ics[k] = datsp[k][-1]
                
				v=[]
				while len(v)<mistakes:
					newEl = np.array([np.random.randint(0,nrow,1),np.random.randint(0,nrow,1)])
					t=add(v,newEl)
					if t:
						v+=[newEl]
				v = np.array(v)
				for k in range(len(v)):
					ics = flip(ics,v[k][0],v[k][1])

				## solve the model
				dat1=run_simulation(tmax,dt,lattice,latticesize,nType,nAmp)
				titleT = "data/mistakes/traj_"+str(nrow)+"x"+str(ncol)+"x1_"+str(nType)+"_n"+str(nAmp)+"_m"+str(mistakes)+"_s"+str(seedInd)+".dat"
				##Put the data into a form that is easier to output
				for k in dat1:
					for i in range(len(dat1[k])):
						dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
				df_Final = pd.DataFrame(dat1)
				df_Final['run'] = seedInd

				df_Final.to_csv(path_or_buf=titleT,index=False)
	    if genDeviations:
		for sdev in sdevList:   
			for seedInd in range(numberOfSims):
				ss = seedList[seedInd]
				np.random.seed(ss)
				nrow=latticesize
				ncol=latticesize
				## Initial condition: Salt and pepper with continuous perturbations (Gaussian added to lattice)
				## First generate the salt and pepper pattern
				datsp=run_simulation(tmax,dt,lattice,latticesize,nType,0,special='SP')

				ics={}
				for k in datsp:
					ics[k] = datsp[k][-1]
					ics[k]=np.random.normal(ics[k],sdev,ics[k].shape)
					ics[k] = (ics[k]>=0)*ics[k]

				## solve the model
				dat1=run_simulation(tmax,dt,lattice,latticesize,nType,nAmp)
				titleT = "data/dev/traj_"+str(nrow)+"x"+str(ncol)+"x1_"+str(nType)+"_n"+str(nAmp)+"_d"+str(sdev)+"_s"+str(seedInd)+".dat"
				##Put the data into a form that is easier to output
				for k in dat1:
					for i in range(len(dat1[k])):
						dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
				df_Final = pd.DataFrame(dat1)
				df_Final['run'] = seedInd

				df_Final.to_csv(path_or_buf=titleT,index=False)
	    if genPatterns:
		    datsp=run_simulation(tmax,dt,lattice,latticesize,nType,0,special='SP')
		    recv ={'N':np.max(datsp['N'][-1]),'D':np.min(datsp['D'][-1]),'I':np.max(datsp['I'][-1])}
		    senv ={'N':np.min(datsp['N'][-1]),'D':np.max(datsp['D'][-1]),'I':np.min(datsp['I'][-1])}
	     	    shapeV = datsp['N'][-1].shape
		    for pattN in pattList:         ## >=0 for pattern
			for seedInd in range(numberOfSims):
				ss = seedList[seedInd]
				np.random.seed(ss)
				nrow=latticesize
				ncol=latticesize
				## Generate the pattern: 
				## 0: only one corner is receivers, 		1: only one corner is senders
				## 2: all senders, 				3: all receivers
				## 4: one line receivers, 			5: one line senders
				## 6: half senders
				## 7: nucleating condition

				## First generate the salt and pepper pattern
				ics=getLatticeICS(pattN,recv,senv,shapeV)

				## solve the model
				dat1=run_simulation(tmax,dt,lattice,latticesize,nType,nAmp)
				titleT = "data/pattern/traj_"+str(nrow)+"x"+str(ncol)+"x1_"+str(nType)+"_n"+str(nAmp)+"_p"+str(pattN)+"_s"+str(seedInd)+".dat"
				##Put the data into a form that is easier to output
				for k in dat1:
					for i in range(len(dat1[k])):
						dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
				df_Final = pd.DataFrame(dat1)
				df_Final['run'] = seedInd

				df_Final.to_csv(path_or_buf=titleT,index=False)

	    if genRandom:
			for seedInd in range(numberOfSims):
				ss = seedList[seedInd]
				np.random.seed(ss)
				nrow=latticesize
				ncol=latticesize

				## solve model for randomized initial conditions
				dat1=run_simulation(tmax,dt,lattice,latticesize,nType,nAmp)
				titleT = "data/random/traj_"+str(nrow)+"x"+str(ncol)+"x1_"+nType+"_n"+str(nAmp)+"_s"+str(seedInd)+".dat"

				##Put the data into a form that is easier to output
				for k in dat1:
					for i in range(len(dat1[k])):
						dat1[k][i] = list(np.around(dat1[k][i].flatten(),decimals=3))
				df_Final = pd.DataFrame(dat1)
				df_Final['run'] = seedInd

				df_Final.to_csv(path_or_buf=titleT,index=False)
