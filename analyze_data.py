## Written by Madeline Galbraith
## Last edited: July 2022

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

from aux_functions import *


## setup to be able to run parallel
from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

## Select all the files we'll be analyzing and 
## separate out to use mpi
## if you don't run parallel all data will be analyzed using a single process
fileList ={}
for i in range(nprocs):
	fileList[i]=[]
count=0
for dirn in os.listdir("data"):
    for filen in os.listdir("data/"+dirn): 
	if ("traj_" in filen):
		fileList[count%nprocs]+=["data/"+dirn+"/"+filen]
		count+=1

### The values in the salt&pepper case
## Sender:    N = 567,    D = 1561, I = 1.2
## Receiver:  N = 5138.8, D = 21.7, I = 802

########### MAIN #############
def main():
	for filen in fileList[myrank]:
			[dN,dD,dI] = get_data(filen,tstart=0)
			[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
                        thresholds={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

			## All data is for each individual time step, lattice, and simulation

			## Calculate and output the similarity of the lattice to checkerboard
			distD =get_results1(dN,dD,thresholds)
			write_Similarity(distD,filen)

			## Calculate and output the avg number of S or R states
    			[avgR,avgS,avglR,avglS,avgER,avgES,states]= get_results2(dN,dD,thresholds)
			write_AvgState(avgR,avgS,avglR,avglS,filen)

			## Calculate and output the correct contacts in the lattice
			[contacts]= get_results3(dN,dD,thresholds)
			write_CorrectContacts(contacts,filen)

			## Calculate and output the time correlation of the lattices
			write_timeCorrelation(dN,dD,thresholds,filen)

			## Calculate and output the transition waiting times for the cells
			[times,neighV,neighS]= get_results7(dN,dD,thresholds,tfin=100000)
			write_WaitingTimes(times,neighV,neighS,filen)

#############################################################
#############################################################
#############################################################
def write_Similarity(res,filen):
    ### write the similarity data to a file  in analysis/folder/simM_[rest the same as data].txt
    title = filen.replace("traj",'simM').replace('data','analysis').replace("dat",'txt')
    fileo = open(title,'w')
    fileo.write("Sim\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s\n" %(res[simX][k]))
    fileo.close()

def write_AvgState(res,res1,res2,res3,filen):
    title = filen.replace("traj",'avgStates').replace('data','analysis').replace("dat",'txt')

    fileo = open(title,'w')
    fileo.write("avgR,avgS,avglR,avglS\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k]))
    fileo.close()

def write_CorrectContacts(res,filen):
    title = filen.replace("traj",'Cont').replace('data','analysis').replace("dat",'txt')
    fileo = open(title,'w')
    fileo.write("contR,contOpp,contS\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX]['SS'])):
	### print it out
	fileo.write("%s,%s,%s\n" %(res[simX]['SS'][k],res[simX]['Opp'][k],res[simX]['RR'][k]))
    fileo.close()
def write_timeCorrelation(ddN,ddD,thresholds,filen):
    states2,states={},{}
    states3={}

    title = filen.replace('traj',"ferr2").replace('data','analysis').replace("dat",'txt')
    fileo = open(title,'w')
    fileo.write("mag,qres\n")
    for sim in ddD:
    	[a,b] = ddN[sim][0].shape
	numberL = a*b
        datD = ddD[sim]
        datN = ddN[sim]        

        threshN={'S':thresholds['S']['N'],'R':thresholds['R']['N']}
        threshD={'S':thresholds['S']['D'],'R':thresholds['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
	compState=states3[sim][0]
	qtmp = compState*states3[sim][:]
	qtmp2 = states3[sim][-1]*states3[sim][:]
	t0=time.time()
        for k in range(len(datD)):
	    if k>0:
	    	qres=(1./numberL)*np.sum(np.mean(qtmp[:k],axis=0))
	    	qres2=(1./numberL)*np.sum(np.mean(qtmp2[:k],axis=0))
	    else:
	    	qres=(1./numberL)*np.sum(np.mean(qtmp[0],axis=0))
	    	qres2=(1./numberL)*np.sum(np.mean(qtmp2[0],axis=0))
	    fileo.write("%s,%s\n" %(qres,qres2))
    fileo.close()

def write_WaitingTimes(time,res1,res2,filen):
    simX = list(time.keys())[0]
    for key in ['effS2R','effR2S']:#,'resS','resR','resIr','resIs']:
        title = filen.replace('traj',"et_eq_"+str(key)).replace('data','analysis').replace("dat",'txt')
        fileo = open(title,'w')
        fileo.write("row,col,time,N,D,R,S\n")
	for rv in time[simX]:
	    for cv in time[simX][rv]:
        	for k in range(len(time[simX][rv][cv][key])):
    	    		fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(rv,cv,time[simX][rv][cv][key][k],res1[simX][rv][cv][key]['N'][k],res1[simX][rv][cv][key]['D'][k],res2[simX][rv][cv][key]['R'][k],res2[simX][rv][cv][key]['S'][k]))
        fileo.close()




#############################
#########################
#########################
#########################
main()
