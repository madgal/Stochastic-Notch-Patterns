import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.stats import pearsonr
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson
from scipy.stats import kde
import matplotlib
matplotlib.rcParams.update({'font.size': 20})


from basicFunctions import *

from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()


### The values in the salt&pepper case
## Sender:    N = 567,    D = 1561, I = 1.2
## Receiver:  N = 5138.8, D = 21.7, I = 802

## The following are thresholds mainly for the 4x4 lattice 
thresholds={'shot':{},'white':{}}

## choose which ones to plot and then plot
## the lattice and similarity together
fileList ={}
for i in range(nprocs):
	fileList[i]=[]


count=0
for filen in os.listdir("data_ND_mcT_trajs_lattices"):
	if ("traj_" in filen):# and (filen in keepRunning):
		fileList[count%nprocs]+=[filen]
		count+=1

########### MAIN #############
def main():

	print(fileList[myrank])
	#fileo = open("contact_info.dat","w")
	#fileo.close()
	for filen in fileList[myrank]:
		if ('traj_' in filen):
			[lattice,noiseType,noiseAmp] = getinfo(filen)
			if True:#noiseType=='shot' and noiseAmp==0:
				print(lattice,noiseType,noiseAmp)
				[dN,dD,dI] = get_data(filen,"data_ND_mcT_trajs_lattices/",tstart=0)
				[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
				if lattice not in thresholds[noiseType].keys():
                                    thresholds[noiseType][lattice]={}
				if noiseAmp not in thresholds[noiseType][lattice].keys():
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				else:
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				writeFullSS1(distD,noiseType,noiseAmp,lattice)

				#[avgR,avgS,avglR,avglS,avgER,avgES,states]= get_results2(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS2(avgR,avgS,avglR,avglS,avgER,avgES,noiseType,noiseAmp,lattice)

				[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				writeFullSS3(contacts,noiseType,noiseAmp,lattice)

def writeFullSS1(res,noiseType,noiseAmp,lattice):
    dir = "figures_lattices/data/"
    simX = list(res.keys())[0]
    title = dir+"simM_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("Sim\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s\n" %(res[simX][k]))
    fileo.close()
def writeFullSS2(res,res1,res2,res3,res4,res5,noiseType,noiseAmp,lattice):
    dir = "figures_lattices/data/"
    simX = list(res.keys())[0]
    title = dir+"avgStates_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("avgR,avgS,avglR,avglS,avgER,avgES\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k],res4[simX][k],res5[simX][k]))
    fileo.close()
def writeFullSS3(res,noiseType,noiseAmp,lattice):
    dir = "figures_lattices/data/"
    simX = list(res.keys())[0]
    title = dir+"Cont_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("contR,contOpp,contS\n")
    for k in range(len(res[simX]['SS'])):
	### print it out
	fileo.write("%s,%s,%s\n" %(res[simX]['SS'][k],res[simX]['Opp'][k],res[simX]['RR'][k]))
    fileo.close()


#########################
#########################
#########################
main()
