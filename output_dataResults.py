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

keepRunning=[ "traj_16x16x1_shot_n13.dat", "traj_16x16x1_white_n130.dat"]


count=0
for filen in os.listdir("data_ND_mcT_trajs_Oct2020"):
	if ("traj_" in filen) and (filen in keepRunning):
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
				[dN,dD,dI] = get_data(filen,"data_ND_mcT_trajs_Oct2020/",tstart=0)
				[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
				if lattice not in thresholds[noiseType].keys():
                                    thresholds[noiseType][lattice]={}
				if noiseAmp not in thresholds[noiseType][lattice].keys():
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				else:
                                    thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
				#distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS1(distD,noiseType,noiseAmp,lattice)

				#[avgR,avgS,avglR,avglS,states]= get_results2(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS2(avgR,avgS,avglR,avglS,noiseType,noiseAmp,lattice)

				#[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS3(contacts,noiseType,noiseAmp,lattice)

				##[mDN,sDN]= get_results4(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				##writeFullSS4(mDN,sDN,noiseType,noiseAmp,lattice)

				#writeFullSS5(dN,dD,thresholds,noiseType,noiseAmp,lattice)
				#writeFullSS5_2(dN,dD,thresholds,noiseType,noiseAmp,lattice)
				#writeFullSS5_f(dN,dD,thresholds,noiseType,noiseAmp,lattice)

				#[switchSR,switchRS,switchiSR,switchiRS]= get_results6(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS6(switchSR,switchRS,switchiSR,switchiRS,noiseType,noiseAmp,lattice)
				#[switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS,times,neighV,neighS]= get_results6(dN,dD,noiseType,noiseAmp,lattice,thresholds)

				#[switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS]= get_results6b(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeFullSS6(switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS,noiseType,noiseAmp,lattice)
				print("Have data")
				[times,neighV,neighS]= get_results7(dN,dD,noiseType,noiseAmp,lattice,thresholds,tfin=100000)
				writeFullSS7(times,neighV,neighS,noiseType,noiseAmp,lattice)

				##print(len(dN))
				##print(dN[0].shape)
				#[simM,avgR,avgS,avglR,avglS,states,contacts,mDN,sDN,mag,qfunc,switchSR,switchRS,switchiSR,switchiRS] = get_results(dN,dD,noiseType,noiseAmp,lattice,thresholds)
				#writeSS(simM,avgR,avgS,avglR,avglS,contacts,mDN,sDN,mag,qfunc,switchSR,switchRS,switchiSR,switchiRS,lattice,noiseType,noiseAmp)
				##[escape_times,neighV,neighS] = getET_states(states,dN,dD)
				##writeETresults(escape_times,neighV,neighS,lattice,noiseType,noiseAmp)

def writeFullSS1(res,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(res.keys())[0]
    title = dir+"simM_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("Sim\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s\n" %(res[simX][k]))
    fileo.close()
def writeFullSS2(res,res1,res2,res3,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(res.keys())[0]
    title = dir+"avgStates_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("avgR,avgS,avglR,avglS\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k]))
    fileo.close()
def writeFullSS3(res,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(res.keys())[0]
    title = dir+"Cont_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("contR,contOpp,contS\n")
    for k in range(len(res[simX]['SS'])):
	### print it out
	fileo.write("%s,%s,%s\n" %(res[simX]['SS'][k],res[simX]['Opp'][k],res[simX]['RR'][k]))
    fileo.close()
def writeFullSS4(res,res1,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(res.keys())[0]
    title = dir+"DN_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("mDN,sDN\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s\n" %(res[simX][k],res1[simX][k]))
    fileo.close()
def writeFullSS5_f(ddN,ddD,thresholds,noiseType,noiseAmp,lattice):
    __thresholds__=thresholds
    states2,states={},{}
    states3={}

    simX = list(ddD.keys())[0]
    dir = "figures_Oct2020/data/"
    title = dir+"ferr_f_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("qresI,qresF\n")
    for sim in ddD:
        datD = ddD[sim]
        datN = ddN[sim]        

        threshN={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['N'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['N']}
        threshD={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['D'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
	compStateI=states3[sim][0]
	compStateF=states3[sim][-1]

        for k in range(len(datD)):
	    qres=np.mean(states3[sim][k]*compStateI)
	    qres1=np.mean(states3[sim][k]*compStateF)
	    fileo.write("%s,%s\n" %(qres,qres1))
    fileo.close()

def writeFullSS5(ddN,ddD,thresholds,noiseType,noiseAmp,lattice):
    __thresholds__=thresholds
    states2,states={},{}
    states3={}

    simX = list(ddD.keys())[0]
    dir = "figures_Oct2020/data/"
    title = dir+"ferr_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("mag,qres\n")
    for sim in ddD:
    	[a,b] = ddN[sim][0].shape
	numberL = a*b
        datD = ddD[sim]
        datN = ddN[sim]        

        threshN={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['N'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['N']}
        threshD={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['D'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
	qtmp = states3[sim][0]*states3[sim][:]
	t0=time.time()
        for k in range(len(datD)):
            mag=getMagnetization(datN[k],datD[k])
	    if k>0:
	    	qres=(1./numberL)*np.sum(np.mean(qtmp[:k],axis=0))
	    else:
	    	qres=(1./numberL)*np.sum(np.mean(qtmp[0],axis=0))
	    fileo.write("%s,%s\n" %(mag,qres))
    fileo.close()
def writeFullSS5_2(ddN,ddD,thresholds,noiseType,noiseAmp,lattice):
    __thresholds__=thresholds
    states2,states={},{}
    states3={}

    simX = list(ddD.keys())[0]
    dir = "figures_Oct2020/data/"
    title = dir+"ferr2_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("mag,qres\n")
    for sim in ddD:
    	[a,b] = ddN[sim][0].shape
	numberL = a*b
        datD = ddD[sim]
        datN = ddN[sim]        

        threshN={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['N'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['N']}
        threshD={'S':__thresholds__[noiseType][lattice][noiseAmp]['S']['D'],'R':__thresholds__[noiseType][lattice][noiseAmp]['R']['D']}
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


def writeFullSS6(res,res1,res2,res3,res4,res5,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(res.keys())[0]
    title = dir+"switches_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
    fileo = open(title,'w')
    fileo.write("switchSR,switchRS,switchiSR,switchiRS,switcheffSR,switcheffRS\n")
    for k in range(len(res[simX])):
	### print it out
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(res[simX][k],res1[simX][k],res2[simX][k],res3[simX][k],res4[simX][k],res5[simX][k]))
    fileo.close()
def writeFullSS7(time,res1,res2,noiseType,noiseAmp,lattice):
    dir = "figures_Oct2020/data/"
    simX = list(time.keys())[0]
    for key in ['effS2R','effR2S','resS','resR','resIr','resIs']:
    	title = dir+"et_eq_"+str(key)+"_"+str(lattice)+"_"+str(noiseType)+"_n"+str(noiseAmp)+"_s"+str(simX)+".txt"
        fileo = open(title,'w')
        fileo.write("row,col,time,N,D,R,S\n")
	for rv in time[simX]:
	    for cv in time[simX][rv]:
        	for k in range(len(time[simX][rv][cv][key])):
    	    		fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(rv,cv,time[simX][rv][cv][key][k],res1[simX][rv][cv][key]['N'][k],res1[simX][rv][cv][key]['D'][k],res2[simX][rv][cv][key]['R'][k],res2[simX][rv][cv][key]['S'][k]))
        fileo.close()



#############################
def writeETresults(et,neighV,neighS,lattice,noiseType,noiseAmp):
    #	trajs[sim][rv][cv]['iR']['D']=DV
    ## Sender:    N = 567,    D = 1561, I = 1.2
    ## Receiver:  N = 5138.8, D = 21.7, I = 802
    ### plot the landscape    
    #times[sim][rv][cv]=
    totals={}
    for sim in et:
        totals[sim]={'effS2R':[],'effR2S':[],'eS2R':[],'eR2S':[],'bS2R':[],'bR2S':[],'etS2R':[],'etR2S':[],'resS':[],'resR':[],'resIr':[],'resIs':[]}

        keylist = totals[sim].keys()
        for key in keylist:
            xval,yval=[],[]
            xva2,yva2=[],[]
            for rv in et[sim]:
                for cv in et[sim][rv]:
                    totals[sim][key]+= et[sim][rv][cv][key]
                    xval+= neighV[sim][rv][cv][key]['D']
                    yval+= neighV[sim][rv][cv][key]['N']
                    xva2+= neighS[sim][rv][cv][key]['S']
                    yva2+= neighS[sim][rv][cv][key]['R']

	    xval= np.array(xval)
	    yval= np.array(yval)
	    xva2= np.array(xva2)
	    yva2= np.array(yva2)
	    totals[sim][key]= np.array(totals[sim][key])
	    ix = np.argwhere(np.isnan(xval)*1. ==1.)[:,0]
	    iy = np.argwhere(np.isnan(yval)*1. ==1.)[:,0]
	    ix2= np.argwhere(np.isnan(xva2)*1. ==1.)[:,0]
	    iy2= np.argwhere(np.isnan(yva2)*1. ==1.)[:,0]
	    ie = np.argwhere(np.isnan(totals[sim][key])*1. ==1.)[:,0]
            ind = np.concatenate([ix,ie,iy])
            ind = np.unique(ind)
            xval = np.delete(xval,ind)
            yval = np.delete(yval,ind)
            xva2 = np.delete(xva2,ind)
            yva2 = np.delete(yva2,ind)
            totals[sim][key] = np.delete(totals[sim][key],ind)

            dir = "figures_Oct2020/data/"
            title = dir+"escapeResults_"+str(lattice)+"_"+str(noiseType)+str(noiseAmp)+"_s"+str(sim)+"_"+str(key)+".txt"
	    fileo = open(title,'w')
	    fileo.write("Next,Dext,Sext,Rext,tescape\n")
	    for i in range(len(totals[sim][key])):
		fileo.write("%s,%s,%s\n" %(xval[i],yval[i],xva2[i],yva2[i],totals[sim][key][i]))
	    fileo.close()

#############################
def writeSS(simM,avgR,avgS,avglR,avglS,contacts,mDN,sDN,mag,qfunc,switchSR,switchRS,switchiSR,switchiRS,lattice,noiseType,noiseAmp):
    ### Plot S/R comp ad Similiarity metric
     for sim in avgR:
	dir = "figures_Oct2020/data/"
        title = dir+"simM_srA_"+str(lattice)+"_"+str(noiseType)+str(noiseAmp)+"_s"+str(sim)+".txt"

	fileo = open(title,'w')
	fileo.write("Sim,avgR,avgS,avglR,avglS,contSS,contOpp,contRR,mDeltaNotch,sDeltaNotch,mag,qres,switchSR,switchRS,switchiSR,switchiRS\n")
	for i in range(len(avgR[sim])):
	    if i==0:
	        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(simM[sim][i],avgR[sim][i],avgS[sim][i],avglR[sim][i],avglS[sim][i],contacts[sim]['SS'][i],contacts[sim]['Opp'][i],contacts[sim]['RR'][i],mDN[sim][i],sDN[sim][i],mag[sim][i],qfunc[sim][i],'','','',''))
	    else:
	        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(simM[sim][i],avgR[sim][i],avgS[sim][i],avglR[sim][i],avglS[sim][i],contacts[sim]['SS'][i],contacts[sim]['Opp'][i],contacts[sim]['RR'][i],mDN[sim][i],sDN[sim][i],mag[sim][i],qfunc[sim][i],switchSR[sim][i-1],switchRS[sim][i-1],switchiSR[sim][i-1],switchiRS[sim][i-1]))
	fileo.close()




#########################
#########################
#########################
main()
