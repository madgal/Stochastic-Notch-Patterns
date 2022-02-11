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
import time


from basicFunctions import *

from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()


### The values in the salt&pepper case
## Sender:    N = 567,    D = 1561, I = 1.2
## Receiver:  N = 5138.8, D = 21.7, I = 802

## The following are thresholds mainly for the 4x4 lattice 
'''
thresholds = {'shot':{0:{'S':{'N':751.14,'D':1190.47},'R':{'N':3764.6,'D':66.95}},
                      4:{'S':{'N':751.14,'D':1190.47},'R':{'N':3764.6,'D':66.95}},
                      8:{'S':{'N':945.63,'D':1061.01},'R':{'N':3355.21,'D':66.95}},
                      16:{'S':{'N':1190.47,'D':596.65},'R':{'N':2116.99,'D':59.66}},
                      32:{'S':{'N':1190.47,'D':596.65},'R':{'N':2116.99,'D':59.66}}
                     },
              'white':{0:{'S':{'N':669.45,'D':1335.73},'R':{'N':3764.6,'D':168.16}},
                      40:{'S':{'N':669.45,'D':1335.73},'R':{'N':3764.6,'D':168.16}},
                      80:{'S':{'N':945.63,'D':1190.47},'R':{'N':2990.33,'D':335.52}},
                      160:{'S':{'N':1190.47,'D':751.14},'R':{'N':1335.73,'D':669.45}},
                      320:{'S':{'N':1190.47,'D':751.14},'R':{'N':1335.73,'D':669.45}}
                      }
             }
'''
thresholds={'shot':{},'white':{}}

## choose which ones to plot and then plot
## the lattice and similarity together
fileList ={}
fileDir ={}
for i in range(nprocs):
	fileList[i]=[]
	fileDir[i]=[]

count=0
'''
for filen in os.listdir("data_ND_mcT_trajs_mar2021/dev"):
	if ("traj_" in filen)  and ("s0" not  in filen) and (('d50_' in filen) or ('d1500_' in filen) or ('d5000_' in filen)) :#
	    	#if (("shot" in filen) and (("n0_" in filen) or ("n5_" in filen) or ("n10_" in filen) or ("n15_" in filen) or ("n20_" in filen))) or (("white" in filen) and ( ("n0" in filen) or ("n50_" in filen) or ("n100_" in filen) or ("n150_" in filen) or ("n200_" in filen))):
		fileList[count%nprocs]+=[filen]
		fileDir[count%nprocs]+=['dev/']
		count+=1
for filen in os.listdir("data_ND_mcT_trajs_mar2021/mistakes"):
	if ("traj_" in filen) and ("s0" not in filen) and (('m13_' in filen) or ('m64_' in filen) or ('m230_' in filen)):# 
		fileList[count%nprocs]+=[filen]
		fileDir[count%nprocs]+=['mistakes/']
		count+=1
'''
for filen in os.listdir("data_ND_mcT_trajs_mar2021/patt"):
	if ("traj_" in filen) and (('p8_' in filen) or ('p7_' in filen)):
		fileList[count%nprocs]+=[filen]
		fileDir[count%nprocs]+=['patt/']
		count+=1

########### MAIN #############
def main():

	#fileo = open("contact_info.dat","w")
	#fileo.close()
	for i in range(len(fileList[myrank])):
		filen = fileList[myrank][i]
		if ('traj_' in filen):
                        [lattice,noiseType,noiseAmp,ext,sim] = getinfoDetN2(filen)
			selected=False#True#
        		#titleO = "figures_mar2021/"+fileDir[myrank][i]+"/data/et_resR_16x16_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_s"+str(sim)+".txt"
			#if (not os.path.exists(titleO)):
			#	selected=True
			#	print filen

 			if True:#selected:#
				#sss = titleO+"\n"+str(lattice)+","+str(ext)+","+str(noiseAmp)+","+str(noiseType)+","+str(sim)
				#print  sss
				print filen

				[dN,dD,dI] = get_data(filen,"data_ND_mcT_trajs_mar2021/"+fileDir[myrank][i],tstart=0)
				try:
					[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
				except:
					[nr,dr,ns,ds,ntr,dtr,nts,dts] = [5138.8,21.7,567.,1561.,2853.,791.,2853.,781.]
				if lattice not in thresholds[noiseType].keys():
                                    thresholds[noiseType][lattice]={}
				if noiseAmp not in thresholds[noiseType][lattice].keys():
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				else:
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				print("Have data")

				#try:
				#	distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#	writeFullSS1(distD,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error sim ",filen

				#try:
				#	[avgR,avgS,avglR,avglS,states]= get_results2(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#	writeFullSS2(avgR,avgS,avglR,avglS,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error sim ",filen

				#try:
				#	[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#	writeFullSS3(contacts,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error sim ",filen
				#try:
				#	[mDN,sDN]= get_results4(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#	writeFullSS4(mDN,sDN,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error DN ",filen

				#try:
				#	[mag,qres]= get_results5(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#	writeFullSS5(mag,qres,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error qf ",filen
				#try:

				[switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS]= get_results6b(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				writeFullSS6(switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)

				try:
					[times,neighV,neighS]= get_results7(dN,dD,noiseType,noiseAmp,lattice,thresholds)
					writeFullSS7(times,neighV,neighS,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				except:
					print "Error sim ",filen
				#try:
				#	filec="traj_16x16x1_"+str(noiseType)+"_n"+str(noiseAmp)+".dat"
				#	[dNc,dDc,dIc] = get_data(filec,"data_ND_mcT_trajs_Oct2020/",tstart=149900)
				#	[nr,dr,ns,ds,ntrc,dtrc,ntsc,dtsc] = get_thresholds(dNc,dDc)
                                #        threshNc={'S':ntsc,'R':ntrc}
                                #        threshDc={'S':dtsc,'R':dtrc}
                		#	dDc = dDc[0]
				#	dNc = dNc[0]        
                                #        [states,states2,states3] = getSimStates(dDc,dNc,threshDc,threshNc)    
				#	comp = states3[-1]
				#	writeFullSS5_2(dN,dD,comp,thresholds,noiseType,noiseAmp,lattice,fileDir[myrank][i],ext,sim)
				#except:
				#	print "Error qf2 ", filen

def writeFullSS5_2(ddN,ddD,compState,thresholds,noiseType,noiseAmp,lattice,Dir,ext,simX):
    __thresholds__=thresholds
    states2,states={},{}
    states3={}

    simX = list(ddD.keys())[0]
    dir = "figures_mar2021/"+Dir+"data/"
    title = dir+"ferr2_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("qres,qres2\n")
    sim = list(ddD.keys())[0]
    [a,b] = ddN[sim][0].shape
    numberL = a*b
    datD = ddD[sim]
    datN = ddN[sim]        

    threshN={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['N'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['N']}
    threshD={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['D'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['D']}
    [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
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
    o.close()


def writeFullSS1(res,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"simM_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("Sim\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s\n" %(res[simX][k]))
    fileo.close()
def writeFullSS2(res,res1,res2,res3,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"avgStates_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("avgR,avgS,avglR,avglS\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k]))
    fileo.close()
def writeFullSS3(res,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"Cont_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("contR,contOpp,contS\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX]['SS'])):
	### print it out
	fileo.write("%s,%s,%s\n" %(res[simX]['SS'][k],res[simX]['Opp'][k],res[simX]['RR'][k]))
    fileo.close()
def writeFullSS4(res,res1,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"DN_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("mDN,sDN\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s\n" %(res[simX][k],res1[simX][k]))
    fileo.close()
def writeFullSS5(res,res1,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"ferr_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("mag,qres\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s\n" %(res[simX][k],res1[simX][k]))
    fileo.close()
def writeFullSS6(res,res1,res2,res3,res4,res5,noiseType,noiseAmp,lattice,Dir,ext,sim):
    dir = "figures_mar2021/"+Dir+"/data/"
    title = dir+"switches_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
    fileo = open(title,'w')
    fileo.write("switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS\n")
    simX = list(res.keys())[0]
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k],res4[simX][k],res5[simX][k]))
    fileo.close()
def writeFullSS7(time,res1,res2,noiseType,noiseAmp,lattice,Dir,ext,sim):
    '''
                            times[sim][rv][cv]['resS']+=[stTime/10.]
                            neighV[sim][rv][cv]['resS']['N']+=[np.mean(tmpNs)]
                            neighV[sim][rv][cv]['resS']['D']+=[np.mean(tmpDs)]
                            neighS[sim][rv][cv]['resS']['R']+=[np.mean(numRs)]
                            neighS[sim][rv][cv]['resS']['S']+=[np.mean(numSs)]
                            effS2R','effR2S','resS','resR','resIr','resIs'
    '''
    dir = "figures_mar2021/"+Dir+"/data/"
    for key in ['effS2R','effR2S','resS','resR','resIr','resIs']:
        title = dir+"et_"+str(key)+"_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_"+str(ext)+"_"+str(sim)+".txt"
        fileo = open(title,'w')
        fileo.write("row,col,time,N,D,R,S\n")
        simX = list(time.keys())[0]
	for rv in time[simX]:
	    for cv in time[simX][rv]:
        	for k in range(len(time[simX][rv][cv][key])):
    	    		fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(rv,cv,time[simX][rv][cv][key][k],res1[simX][rv][cv][key]['N'][k],res1[simX][rv][cv][key]['D'][k],res2[simX][rv][cv][key]['R'][k],res2[simX][rv][cv][key]['S'][k]))
        fileo.close()

#########################
#########################
#########################
main()
