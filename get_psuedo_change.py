import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib
from basicFunctions import *
from plot_script_functions import get_dataVersion2
from plot_script_functions import pseudo_potential
from plot_script_functions import find_minima


def find_thresholds(x, U, i_s, j_s, i_r, j_r,scale):
    #written by Federico Bocci
    '''
    find distances from S and R minima where the probability is decreased by a 10-fold
    thus the pseudopotential increases by a unit
    '''
    ### find __thresholds__:
    # from Receiver: notch moves left, delta moves up
    Uref = U[i_s][j_s]
    i = i_s - 1
    while (U[i][j_s] - Uref) < np.log10(10.*scale):#
        i = i - 1
    notch_thr_R = x[i]
    j = j_s + 1
    while (U[i_s][j] - Uref) <  np.log10(10.*scale):
        j = j + 1
    delta_thr_R = x[j]
    # from Sender: notch moves right, delta moves down
    Uref = U[i_r][j_r]
    i = i_r + 1
    while (U[i][j_r] - Uref) <  np.log10(10.*scale):
        i = i + 1
    notch_thr_S = x[i]
    j = j_r - 1
    while (U[i_r][j] - Uref) <  np.log10(10.*scale):
        j = j - 1
    delta_thr_S = x[j]
    return notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S
def get_thresholds(dN,dD,scale):
    ## based off code written by Federico for finding the minima and getting the threshold at 0.1 of that.

    N= dN#[0]
    D= dD#[0]

    bins = np.logspace(np.log10(0.1),np.log10(100000),num=121,base=10)
    x = np.zeros(bins.size-1)
    for i in range(x.size):
        x[i] = (bins[i]+bins[i+1])/2.

    U = pseudo_potential(N,D,bins)
    i_s,j_s,i_r,j_r = find_minima(U,bins,x,x)
    ntr,dtr,nts,dts = find_thresholds(x,U,i_s,j_s,i_r,j_r,scale)
    nr = bins[i_s]
    dr = bins[j_s]
    ns = bins[i_r]
    ds = bins[j_r]

    #print('nr,dr,ns,ds')
    #print(nr,dr,ns,ds)
    #print('ntr,dtr,nts,dts')
    #print(ntr,dtr,nts,dts)

    return [nr,dr,ns,ds,ntr,dtr,nts,dts]


##################3333
##################3333
def get_pseudo(nt,na,scale):
    # parameters for plotting
    sim = 20
    #sim_points = 10000
    sim_points = 20000#150000
    #n_eq = 1000
    n_eq = 10000
    bins = np.logspace(np.log10(0.1), np.log10(100000), num=61, base=10)
    ncell = 16

    N, D = get_dataVersion2('data/random/traj_16x16x1_'+str(nt)+'_n'+str(na)+'.dat',10000,25000,0)
    U = pseudo_potential(N, D, bins)
    
    x = np.zeros(bins.size - 1)
    for i in range(x.size):
        x[i] = (bins[i] + bins[i + 1]) / 2.
        
    notch=bins[0:-1]
    delta=bins[0:-1]
    xint=[10, 30000]
    yint=[0.1, 10000]

    # find minima
    i_s, j_s, i_r, j_r = find_minima(U, bins, x, x)
    notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S = find_thresholds(x, U, i_s, j_s, i_r, j_r,scale)

    
    return [notch,delta,np.transpose(U),xint,yint,[bins[i_s], bins[i_r]], [bins[j_s], bins[j_r]],
	     notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S]	

#### Figure 5
'''
[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',0,0.9)
fileo  = open("main_data/figPP_shot0_lower.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()


[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',10,0.9)
fileo  = open("main_data/figPP_shot10_lower.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',0,1.1)
fileo  = open("main_data/figPP_shot0_higher.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()


[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',10,1.1)
fileo  = open("main_data/figPP_shot10_higher.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()
####
[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',0,0.9)
fileo  = open("main_data/figPP_white0_lower.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()


[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',100,0.9)
fileo  = open("main_data/figPP_white1000_lower.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',0,1.1)
fileo  = open("main_data/figPP_white0_higher.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()


[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',100,1.1)
fileo  = open("main_data/figPP_white1000_higher.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima,ntr,dtr,nts,dts\n")
count=0
for i in range(len(N)):
        for j in range(len(D)):
                if count==0:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0],ntr,dtr,nts,dts))
                elif count==1:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1],-1,-1,-1,-1))
                else:
                        fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1,-1,-1,-1,-1))
                count+=1
fileo.close()
'''


########################
########################
########################
simX=0
thresholds={'shot':{},'white':{}}
dirN='data/random/'
filen = dirN+'traj_16x16x1_shot_n0.dat'
dN, dD = get_dataVersion2(filen,0,-1,0)
#[dN,dD,dI] = get_data(filen,tstart=0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,0.9)
print "shot 0  t=0.9\t", nr,dr,ns,ds,ntr,dtr,nts,dts

lattice='16x16'
noiseType='shot'
noiseAmp=0
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}


[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_shot0_lower.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()
fileo = open("main_data/dat_cont_shot0_lower.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,1.1)
print "shot 0  t=1.1\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_shot0_higher.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_shot0_higher.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()


thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,2.)
print "shot 0  t=2\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,dirN,tstart=0)
fileo = open("main_data/dat_sim_shot0_high20.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_shot0_high20.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

########################
dirN='data/random/'
filen = dirN+'traj_16x16x1_shot_n10.dat'
thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,0.9)
print "shot 10  t=0.9\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_shot10_lower.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_shot10_lower.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,1.1)
print "shot 10  t=1.1\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_shot10_higher.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_shot10_higher.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()


thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,2.)
print "shot 10  t=2\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_shot10_high20.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_shot10_high20.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

########################
dirN='data/random/'
filen =dirN+ 'traj_16x16x1_white_n100.dat'
thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,0.9)
print "white 100  t=0.9\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_white1000_lower.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_white1000_lower.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(dirN+filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,1.1)
print "white 100  t=1.1\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_white1000_higher.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_white1000_higher.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()


thresholds={'shot':{},'white':{}}
dN, dD = get_dataVersion2(filen,0,-1,0)
[nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD,2.)
print "white 100  t=2\t", nr,dr,ns,ds,ntr,dtr,nts,dts
if lattice not in thresholds[noiseType].keys():
       thresholds[noiseType][lattice]={}
if noiseAmp not in thresholds[noiseType][lattice].keys():
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}
else:
       thresholds[noiseType][lattice][noiseAmp]={'S':{'N':nts,'D':dts},'R':{'N':ntr,'D':dtr}}

[dN,dD,dI] = get_data(filen,tstart=0)
fileo = open("main_data/dat_sim_white1000_high20.txt",'w')
distD =get_results1(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(distD[simX])):
	fileo.write("%s\n" %(distD[simX][k]))
fileo.close()

fileo = open("main_data/dat_cont_white1000_high20.txt",'w')
[contacts]= get_results3(dN,dD,noiseType,noiseAmp,lattice,thresholds)
for k in range(len(contacts[0]['SS'])):
	fileo.write("%s\n" %(contacts[simX]['Opp'][k]))
fileo.close()

