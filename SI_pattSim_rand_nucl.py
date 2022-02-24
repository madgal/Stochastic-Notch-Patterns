import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from auxFunctions import *

##Written by Madeline Galbraith

tend=3000

def getrandPatt():
    # returns the 2d pattern as a function of delta in the cells
    # starting from randomized initial conditions
    [dN,dD,dI] = get_data("data/random/traj_16x16x1_white_n0.dat",tstart=0)
    kk = dN.keys()[0]
    [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    
    #iD = dtr+ (dts-dtr)/2.
    #iN = nts+ (ntr-nts)/2.
    #patt=[]
    minv=3000
    maxv=0
    for i in range(tend):#len(dD[kk])):
	if minv>np.min(dD[kk][i]):
		minv = np.min(dD[kk][i])
	if maxv<np.max(dD[kk][i]):
		maxv = np.max(dD[kk][i])
    #	    patt += [(dD[kk][i]>iD)*(dN[kk][i]<iN)*1.+(dD[kk][i]<iD)*(dN[kk][i]>iN)*2.]
    
    return dD[kk][:tend],minv,maxv#patt
def getcheckPatt():
    # returns the 2d pattern as a function of delta in the cells
    # starting from checkerboard initial conditions
    [dN,dD,dI] = get_data("data/pattern/traj_16x16x1_shot_n0_p8_s0.dat",tstart=0)
    kk = dN.keys()[0]
    [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    
    #iD = dtr+ (dts-dtr)/2.
    #iN = nts+ (ntr-nts)/2.
    #patt=[]
    #for i in range(tend):#len(dD[kk])):
    #	    patt += [(dD[kk][i]>iD)*(dN[kk][i]<iN)*1.+(dD[kk][i]<iD)*(dN[kk][i]>iN)*2.]
    
    return dD[kk][:tend],minv,maxv#patt

def getnuclPatt():
    # returns the 2d pattern as a function of delta in the cells
    # starting from nucleating initial conditions
    [dN,dD,dI] = get_data("data/pattern/traj_16x16x1_shot_n0_p7_s0.dat",tstart=0)
    kk = dN.keys()[0]
    [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    
    #iD = dtr+ (dts-dtr)/2.
    #iN = nts+ (ntr-nts)/2.
    #patt=[]
    minv=3000
    maxv=0
    for i in range(tend):#len(dD[kk])):
	if minv>np.min(dD[kk][i]):
		minv = np.min(dD[kk][i])
	if maxv<np.max(dD[kk][i]):
		maxv = np.max(dD[kk][i])
    #        patt += [(dD[kk][i]>iD)*(dN[kk][i]<iN)*(1.)+(dD[kk][i]<iD)*(dN[kk][i]>iN)*(2.)]
    
    return dD[kk][:tend],minv,maxv#patt

def getSim_randFull_det():
    # returns the similarity metric as a function of time
    # starting from randomized initial conditions
    fileN = "data/random/simM_16x16_shot_n0_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    x =np.arange(0,len(filX['Sim'].values[:tend]))*0.1
    return x,filX['Sim'].values[:tend]

def getSim_checkFull_det():
    # returns the similarity metric as a function of time
    # starting from checkerboard initial conditions
    fileN = "data/pattern/simM_16x16_shot_n0_p8_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    x =np.arange(0,len(filX['Sim'].values[:tend]))*0.1
    return x,filX['Sim'].values[:tend]
def getSim_nuclFull_det():
    # returns the similarity metric as a function of time
    # starting from nucleating initial conditions
    fileN = "data/pattern/simM_16x16_shot_n0_p7_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    x =np.arange(0,len(filX['Sim'].values[:tend]))*0.1
    return x,filX['Sim'].values[:tend]



######################
######################
######################
######################
pattIn,minv,maxv=getnuclPatt()
timeN,sim_nucl = getSim_nuclFull_det()

## output all the figures to generate a movie showing the patterns and similarity
## startin from nucleating initial conditions
vmin,vmax = minv,maxv
for k in range(len(timeN)):
	fig2 = plt.figure(constrained_layout=True)
	spec2 = gridspec.GridSpec(ncols=15, nrows=1, figure=fig2)
	axb = fig2.add_subplot(spec2[0, 0:9])
	axd = fig2.add_subplot(spec2[0, 9:14])
	axe = fig2.add_subplot(spec2[0,14])

	axb.set_xlim(0,np.max(timeN))
	axb.set_ylim(0.3,1.05)

	if k<=0:
		axb.plot(timeN[k],sim_nucl[k],'-',linewidth=4)
	else:
		axb.plot(timeN[:k],sim_nucl[:k],'-',linewidth=4)


	im = axd.imshow(pattIn[k],cmap='binary',vmin=vmin,vmax=vmax)#,norm=norm,cmap=cm)#'binary')

	cbar = plt.colorbar(im,cax=axe,label='Delta (molecules)')
	fig2.savefig('movie_si_patt/si_pattnuc_'+str(k)+'.png',bbox_inches='tight')
	#plt.show()
	plt.close()


pattIr,minv,maxv=getrandPatt()
time,sim_rand = getSim_randFull_det()

## output all the figures to generate a movie showing the patterns and similarity
## startin from randomized initial conditions
vmin,vmax = minv,maxv
for k in range(len(timeN)):
	fig2 = plt.figure(constrained_layout=True,figsize=(8,16))
	spec2 = gridspec.GridSpec(ncols=15, nrows=1, figure=fig2)
	axa = fig2.add_subplot(spec2[0, 0:9])
	axc = fig2.add_subplot(spec2[0, 9:14])
	axe = fig2.add_subplot(spec2[0,14])

	axa.set_xlim(0,np.max(timeN))
	axa.set_ylim(0.3,1.05)

	if k<=0:
		axa.plot(time[k],sim_rand[k],'-k',linewidth=4)
	else:
		axa.plot(time[:k],sim_rand[:k],'-k',linewidth=4)

	im =axc.imshow(pattIr[k],cmap='binary',vmin=vmin,vmax=vmax)#,norm=norm,cmap=cm)#'binary')

	cbar = plt.colorbar(im,cax=axe,label='Delta (molecules)')
	fig2.savefig('movie_si_patt/si_pattran_'+str(k)+'.png',bbox_inches='tight')
	plt.close()
