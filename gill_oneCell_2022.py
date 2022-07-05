# One cell Gillespie simulation to quantify intrinsic noise
# Written by Madeline Galbraith
# Last modified July 2022

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json


####################
## Allow the calculation to be parallelized
####################
from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

####################
### The values in the checkerboard
####################
## Sender:    N = 567,    D = 1561, I = 1.2
## Receiver:  N = 5138.8, D = 21.7, I = 802

####################
#### parameters ####
####################
N0 = 5.0e+2      
D0 = 1.0e+3  
I0 = 2.0e+2      ## hill factor

kc = 5.0e-4     #1/hour/molecule
kt = 5.0e-5     #1/hour/molecule

nI = 2.0e+0

lambdaN = 2.0e+0   # no dimension
lambdaD = 0.0e+0   # no dimension

gamma = 1.0e-1   #1/hour
gamma_i = 5.0e-1 #1/hour

####################
####################
#################
def shiftedHill(X,X0,nX,lambdY):
	return lambdY+(1.0 - lambdY)/(1.0+(X/X0)**nX)
####################
def simulate_gill(iterations,ss,Next,Dext,equil=30000):

	## set up the simulation
	np.random.seed(ss)
	time=0.
	results={}
	
	N,D,I = np.random.randint(0,6000,1)[0],np.random.randint(0,2000,1)[0],np.random.randint(0,1000,1)[0]
	results['time']=[time]
	results['N']=[N]
	results['D']=[D]
	results['I']=[I]

	## run gillespie algorithm
	for ind in range(iterations):
		tmp=[]
		r1,r2 = np.random.uniform(0,1,2)
		itt = 1./__a0(N,D,I,Next,Dext)*np.log(1./r1)

	
		[N,D,I] = __updateRes(N,D,I,r2,Next,Dext)
		time+=itt
		results['time']+=[time]
		results['N']+=[N]
		results['D']+=[D]
		results['I']+=[I]


	results['time'] = results['time'][equil:]
	results['N'] = results['N'][equil:]
	results['D'] = results['D'][equil:]
	results['I'] = results['I'][equil:]
	
	return results


#### Functions and definitions directly related to Gillespie algorithm
def __updateRes(N,D,I,r2,Next,Dext):
	if (__X1(N,D,I,Next,Dext)<=r2) and (r2< __X2(N,D,I,Next,Dext)):
		return [N-1,D,I+1]
	elif (__X2(N,D,I,Next,Dext)<=r2) and (r2< __X3(N,D,I,Next,Dext)):
		return [N,D-1,I]
	elif (__X3(N,D,I,Next,Dext)<=r2) and (r2< __X4(N,D,I,Next,Dext)):
		return [N-1,D-1,I]
	elif (__X4(N,D,I,Next,Dext)<=r2) and (r2< __X5(N,D,I,Next,Dext)):
		return [N-1,D,I]
	elif (__X5(N,D,I,Next,Dext)<=r2) and (r2< __X6(N,D,I,Next,Dext)):
		return [N,D-1,I]
	elif (__X6(N,D,I,Next,Dext)<=r2) and (r2< __X7(N,D,I,Next,Dext)):
		return [N,D,I-1]
	elif (__X7(N,D,I,Next,Dext)<=r2) and (r2< __X8(N,D,I,Next,Dext)):
		return [N+1,D,I]
		#return [N,D+1,I]
	else:#if (__X8(N,D,I,Next,Dext)<=r2) and (r2< __X9(N,D,I,Next,Dext)):
		return [N,D+1,I]
		#return [N+1,D,I]


##############
### propensities
##############
def __a1(N,D,I,Next,Dext):
	return kt*N*Dext				#	0.00005*N*Dext*(I-1)  gives N=[0]
def __a2(N,D,I,Next,Dext):
	return kt*D*Next 				#	0.00005*D*Next gives N=[0]
def __a3(N,D,I,Next,Dext):
	return kc*N*D 					#	0.0005*N*D gives N=[0,6000]
def __a4(N,D,I,Next,Dext):
	return gamma*N 					#	0.1*N  gives N=[0,600]
def __a5(N,D,I,Next,Dext):
	return gamma*D					#	0.1*D gives N=[0,200]
def __a6(N,D,I,Next,Dext):
	return gamma_i*I 				#	0.5*I gives N=[0,500]
def __a7(N,D,I,Next,Dext):
	return N0*shiftedHill(I,I0,nI,lambdaN) 		#	500 gives N=[500]
def __a8(N,D,I,Next,Dext):
	return D0*shiftedHill(I,I0,nI,lambdaD) 		#	1000 gives N=[1000]
def __a0(N,D,I,Next,Dext):
	return __a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext)+__a4(N,D,I,Next,Dext)+__a5(N,D,I,Next,Dext)+__a6(N,D,I,Next,Dext)+__a7(N,D,I,Next,Dext)+__a8(N,D,I,Next,Dext)


##############
## probabilities
##############
def __X1(N,D,I,Next,Dext):
	return 0
def __X2(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X3(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X4(N,D,I,Next,Dext):
	return  (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X5(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext)+__a4(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X6(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext)+__a4(N,D,I,Next,Dext)+__a5(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X7(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext)+__a4(N,D,I,Next,Dext)+__a5(N,D,I,Next,Dext)+__a6(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X8(N,D,I,Next,Dext):
	return (__a1(N,D,I,Next,Dext)+__a2(N,D,I,Next,Dext)+__a3(N,D,I,Next,Dext)+__a4(N,D,I,Next,Dext)+__a5(N,D,I,Next,Dext)+__a6(N,D,I,Next,Dext)+__a7(N,D,I,Next,Dext))/__a0(N,D,I,Next,Dext)
def __X9(N,D,I,Next,Dext):
	return 1




###############

##############################
############ Start simulation ###########
##############################
def main():
	seedList=pd.read_csv(filepath_or_buffer='seedsForSims.txt',header=None).values[:,0]
	simNum=10
	iterations=200000
	Next = 567
	Dext = 1561

	simList=[myrank,myrank+nprocs]
	etime=0
	pairs = [[567,1561],[5139,22]]
	for [Next,Dext] in pairs:
		for seedInd in simList:
			ss = seedList[seedInd]

			res = simulate_gill(iterations,ss,Next,Dext,equil=etime)

			fullRes={'res':res,'Next':Next,'Dext':Dext,'equil':etime,'seed':ss,'iterations':iterations}
			title='data_gil/oneCell/data_'+str(Next)+"_"+str(Dext)+"_s"+str(seedInd)
                        json_data=json.dumps(fullRes,indent=4)
                        with open(title+'.json','w') as outfile:
                                outfile.write(json_data)

main()
