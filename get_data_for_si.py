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

def sim_OG(X):
        [nr,nc] = X.shape

        totalC,count=0,0
        for i in range(nr):
                corr = corrF(X[i],X[(i+1)%nr])

                corr = (1.-corr)/2 # convert from (-1,1) to (0,1)
                totalC+=corr
                count+=1.
        totalC=totalC/count

        totalR,count=0,0
        for i in range(nr):
                corr = corrF(X[:,i],X[:,(i+1)%nr])

                corr = (1.-corr)/2 # convert from (-1,1) to (0,1)
                totalR+=corr
                count+=1.
        totalR=totalR/count

        return (totalR+totalC)/2.

def corrF(x,y):
    xbar = np.mean(x)
    ybar = np.mean(y)
    sx=np.std(x)
    sy=np.std(y)    
    if sx==0 or sy==0 or len(x)==0:
        return 1.
    return np.sum((x-xbar)*(y-ybar))/(sx*sy)/len(x)

##############################
##############################
##############################
##############################

def getSpecialPatt(ext):
    [dN,dD,dI] = get_data("data/pattern/traj_16x16x1_"+str(ext)+".dat",tstart=0)
    kk = dN.keys()[0]
    try:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    except:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = [5138.8,21.7,567.,1561.,2853.,791.,2853.,781.]

   	 
    #iD = dtr+ (dts-dtr)/2.
    #iN = nts+ (ntr-nts)/2.
    #pattI = (dD[kk][0]>iD)*(dN[kk][0]<iN)*1.+(dD[kk][0]<iD)*(dN[kk][0]>iN)*2.
    #pattF = (dD[kk][-1]>iD)*(dN[kk][-1]<iN)*1.+(dD[kk][-1]<iD)*(dN[kk][-1]>iN)*2.
    pattI = dD[kk][0]
    pattF = dD[kk][-1]

    datD = dD[dD.keys()[0]]
    datN = dN[dN.keys()[0]]
    distD=[]
    
    threshN={'S':nts,'R':ntr}
    threshD={'S':dts,'R':dtr}
    states,states2,states3 = getSimStates(datD,datN,threshD,threshN)

    cc=[]
    for k in range(len(datD)):
        distD+=[sim_OG(datD[k])]
        tmp = getCont(states2[k])
        cc+=[tmp[1]/5.12]
        
    time = np.arange(0,len(datD))*0.1
    
    return pattI,pattF,distD,cc,time


################################
################################
################################
def getMistakesPatt(ext,nt,na):
    
    if na==0:
        [dN,dD,dI] = get_data("data/mistakes/traj_16x16x1_"+str(ext)+"_s0.dat",tstart=0)
    else:
        [dN,dD,dI] = get_data("data/mistakes/traj_16x16x1_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.dat",tstart=0)
    kk = dN.keys()[0]

    try:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    except:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = [5138.8,21.7,567.,1561.,2853.,791.,2853.,781.]

    
    iD = dtr+ (dts-dtr)/2.
    iN = nts+ (ntr-nts)/2.
    pattI = dD[kk][0]#(dD[kk][0]>iD)*(dN[kk][0]<iN)*1.+(dD[kk][0]<iD)*(dN[kk][0]>iN)*2.
    pattF = dD[kk][-1]#(dD[kk][-1]>iD)*(dN[kk][-1]<iN)*1.+(dD[kk][-1]<iD)*(dN[kk][-1]>iN)*2.

    datD = dD[dD.keys()[0]]
    datN = dN[dN.keys()[0]]
    distD=[]
    
    threshN={'S':nts,'R':ntr}
    threshD={'S':dts,'R':dtr}
    states,states2,states3 = getSimStates(datD,datN,threshD,threshN)

    cc=[]
    for k in range(len(datD)):
        distD+=[sim_OG(datD[k])]
        tmp = getCont(states2[k])
        cc+=[tmp[1]/5.12]
        
    time = np.arange(0,len(datD))*0.1
    
    return pattI,pattF,distD,cc,time


#######################
#######################
#######################
def getDevPatt(ext,nt,na):
    
    if na==0:
        [dN,dD,dI] = get_data("data/dev/traj_16x16x1_"+str(ext)+"_s0.dat",tstart=0)
    else:
        [dN,dD,dI] = get_data("data/dev/traj_16x16x1_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.dat",tstart=0)
    kk = dN.keys()[0]

    try:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)
    except:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = [5138.8,21.7,567.,1561.,2853.,791.,2853.,781.]

    
    iD = dtr+ (dts-dtr)/2.
    iN = nts+ (ntr-nts)/2.
    pattI =dD[kk][0]# (dD[kk][0]>iD)*(dN[kk][0]<iN)*1.+(dD[kk][0]<iD)*(dN[kk][0]>iN)*2.
    pattF = dD[kk][-1]#(dD[kk][-1]>iD)*(dN[kk][-1]<iN)*1.+(dD[kk][-1]<iD)*(dN[kk][-1]>iN)*2.

    datD = dD[dD.keys()[0]]
    datN = dN[dN.keys()[0]]
    distD=[]
    
    threshN={'S':nts,'R':ntr}
    threshD={'S':dts,'R':dtr}
    states,states2,states3 = getSimStates(datD,datN,threshD,threshN)

    cc=[]
    for k in range(len(datD)):
        distD+=[sim_OG(datD[k])]
        tmp = getCont(states2[k])
        cc+=[tmp[1]/5.12]
    time = np.arange(0,len(datD))*0.1
    
    return pattI,pattF,distD,cc,time


##################3333
##################3333
def get_pseudo(nt,na):
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
    notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S = find_thresholds(x, U, i_s, j_s, i_r, j_r)

    
    return [notch,delta,np.transpose(U),xint,yint,[bins[i_s], bins[i_r]], [bins[j_s], bins[j_r]],
	     notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S]	

##################3333
##################3333
##################3333
##################3333

def getAvgs_rand():
    tstart=10000
    tfin=100000
    xa,xb=[],[]
    ya1,ya2,yb1,yb2=[],[],[],[]
    yera1,yera2,yerb1,yerb2=[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp,tmp1,tmp2,tmp3=[],[],[],[]
            for sim in range(20):
                try:
                    fileN = "analysis/random/avgStates_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['avgER'].values[tstart:tfin])
                    tmp1+=list(filX['avgES'].values[tstart:tfin])
                except:
                    count=0
	    tmp = np.array(tmp)
	    tmp1 = np.array(tmp1)
            if nt=='white':
                xa+=[na*10.]
                ya1+=[np.mean(tmp)]
                ya2+=[np.mean(tmp1)]
                yera1+=[np.std(tmp)]
                yera2+=[np.std(tmp1)]
            elif nt=='shot':
                xb+=[na]
                yb1+=[np.mean(tmp)]
                yb2+=[np.mean(tmp1)]
                yerb1+=[np.std(tmp)]
                yerb2+=[np.std(tmp1)]
    return xa,ya1,yera1,xb,yb1,yerb1,ya2,yb2,yera2,yerb2

############################33
############################33
############################33
def getAvgs_check():
    tstart=10000
    xa,xb=[],[]
    ya1,ya2,ya3,ya4,yb1,yb2,yb3,yb4=[],[],[],[],[],[],[],[]
    yera1,yera2,yera3,yera4,yerb1,yerb2,yerb3,yerb4=[],[],[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp,tmp1,tmp2,tmp3=[],[],[],[]
            for sim in range(20):
                try:
                    fileN = "analysis/pattern/avgStates_16x16_"+nt+"_n"+str(na)+"_p8_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['avgER'].values[tstart:])
                    tmp1+=list(filX['avgES'].values[tstart:])
                except:
                    count=0
            if nt=='white':
                xa+=[na*10.]
                ya1+=[np.mean(tmp)]
                ya2+=[np.mean(tmp1)]
                yera1+=[np.std(tmp)]
                yera2+=[np.std(tmp1)]
            elif nt=='shot':
                xb+=[na]
                yb1+=[np.mean(tmp)]
                yb2+=[np.mean(tmp1)]
                yerb1+=[np.std(tmp)]
                yerb2+=[np.std(tmp1)]
    return xa,ya1,yera1,xb,yb1,yerb1,ya2,yb2,yera2,yerb2

#############
##################################3
##################################3
def getCont_rand_ind(sim):
    tstart=10000
    tfin=100000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            try:
                    fileN = "analysis/random/Cont_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN)#.dropna()
                    ind= np.argwhere(np.isnan(filX['contOpp']))[:,0]
                    filX['contOpp'][ind]=0
                    tmp+=list(filX['contOpp'].values[tstart:tfin]/5.12)
            except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
    return xa,ya,yera,xb,yb,yerb

def getSim_rand_ind(sim):
    tstart=10000
    tfin=100000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            try:
                    fileN = "analysis/random/simM_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN)#.dropna()
                    ind= np.argwhere(np.isnan(filX['Sim']))[:,0]
                    filX['Sim'][ind]=0
                    tmp+=list(filX['Sim'].values[tstart:tfin])
            except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
    return xa,ya,yera,xb,yb,yerb

def getSim_rand():
    tstart=10000
    tfin=100000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "analysis/random/simM_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN)#.dropna()
                    ind= np.argwhere(np.isnan(filX['Sim']))[:,0]
                    filX['Sim'][ind]=0
                    tmp+=list(filX['Sim'].values[tstart:tfin])
                except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
    return xa,ya,yera,xb,yb,yerb
def getSimSP():
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "analysis/pattern/simM_16x16_"+nt+"_n"+str(na)+"_p8_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN)#.dropna()
                    ind= np.argwhere(np.isnan(filX['Sim']))[:,0]
                    filX['Sim'][ind]=0
                    tmp+=list(filX['Sim'].values[tstart:])
                except:
                    count=0
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
                
    return xa,ya,yera,xb,yb,yerb
##################3
def getContacts(Dir,ext):
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "analysis/"+Dir+"/Cont_16x16_"+nt+"_n"+str(na)+"_"+ext+"_s0.txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['contOpp'].values[tstart:]/5.12)
                except:
                    count=0
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
                
    return xa,ya,yera,xb,yb,yerb

def getSim(Dir,ext):
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "analysis/"+Dir+"/simM_16x16_"+nt+"_n"+str(na)+"_"+ext+"_s0.txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['Sim'].values[tstart:])
                except:
                    count=0
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
                
    return xa,ya,yera,xb,yb,yerb

##################3
##################3
##################3
##################3
def getSim_lattice(ext):
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "data/lattices/simM_"+ext+"_"+nt+"_n"+str(na)+"_s0.txt"
                    filX=pd.read_csv(fileN).dropna() 
                    tmp+=list(filX['Sim'].values[tstart:])
                except:
                    count=0
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp)]
                yera+=[np.std(tmp)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp)]
                yerb+=[np.std(tmp)]
                
    return xa,ya,yera,xb,yb,yerb

def getContacts_lattice(ext,scale):
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "data/lattices/Cont_"+ext+"_"+nt+"_n"+str(na)+"_s0.txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['contOpp'].values[tstart:])
                except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp/scale)]
                yera+=[np.std(tmp/scale)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp/scale)]
                yerb+=[np.std(tmp/scale)]
                
    return xa,ya,yera,xb,yb,yerb


############
############
############
def getET_avg(df):
    x,y = np.round(df['D'],0),np.round(df['N'],-1)
    yval=(np.unique(y))
    xval=(np.unique(x))
    #plt.plot(yval)
    #plt.plot(xval)

    print len(xval),len(yval)
    a,b,c=[],[],[]
    std=[]
    bk =False

    count=0
    for el in xval:
        inds = np.argwhere(x==el)[:,0]
        for el2 in yval:
            ind2 = np.argwhere(y==el2)[:,0]
            ind = np.intersect1d(inds,ind2)
            if len(ind)>0:
                a+=[el]
                b+=[el2]
                c+=[np.mean(df['time'].values[ind])]
                std+=[np.std(df['time'].values[ind])]

            count+=1



    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    std = np.array(std)
    return [a,b,c,std]

###################
###################
###################
###################
def getETvSR_avg(df):
    x,y = np.round(df['S'],3),np.round(df['R'],3)
    yval=(np.unique(y))
    xval=(np.unique(x))
    #plt.plot(yval)
    #plt.plot(xval)

    print len(xval),len(yval)
    a,b,c=[],[],[]
    std=[]
    bk =False

    count=0
    for el in xval:
        inds = np.argwhere(x==el)[:,0]
        for el2 in yval:
            ind2 = np.argwhere(y==el2)[:,0]
            ind = np.intersect1d(inds,ind2)
            if len(ind)>0:
                a+=[el]
                b+=[el2]
                c+=[np.mean(df['time'].values[ind])]
                std+=[np.std(df['time'].values[ind])]

            count+=1

    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    std = np.array(std)
    return [a,b,c,std]

#####################
#####################
#####################
#####################
#####################
def getSim_checkFull_noise():
    tfin=10000
    fileN = "analysis/pattern/simM_16x16_white_n0_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n0w=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_white_n50_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n50=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_white_n130_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n130=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_white_n200_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n200=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_shot_n0_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n0s=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_shot_n5_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n5=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_shot_n13_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n13=filX['Sim'][:tfin]
    fileN = "analysis/pattern/simM_16x16_shot_n20_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n20=filX['Sim'][:tfin]

    x =np.arange(0,len(filX['Sim'][:tfin]))*0.1
    return x,n0w,n0s,n50,n130,n200,n5,n13,n20

def getCont_checkFull_noise():
    tfin=10000
    fileN = "analysis/pattern/Cont_16x16_white_n0_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n0w=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_white_n50_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n50=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_white_n130_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n130=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_white_n200_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n200=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_shot_n0_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n0s=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_shot_n5_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n5=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_shot_n13_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n13=filX['contOpp'][:tfin]/5.12
    fileN = "analysis/pattern/Cont_16x16_shot_n20_p8_s0.txt"
    filX=pd.read_csv(fileN)
    n20=filX['contOpp'][:tfin]/5.12

    x =np.arange(0,len(filX['contOpp'][:tfin]))*0.1
    return x,n0w,n0s,n50,n130,n200,n5,n13,n20

####################33
####################33
##################33
def qcorr_rand(nt,na):
    tfin=10000
    fileN= "analysis/random/ferr2_16x16_"+nt+"_n"+str(na)+"_s0.txt"
    fileN2= "analysis/random/ferr_16x16_"+nt+"_n"+str(na)+"_s0.txt"
    filX = pd.read_csv(fileN).dropna()
    fil2 = pd.read_csv(fileN2).dropna()

    yv1=fil2['qres'][0:tfin]
    xv1 = np.arange(0,len(yv1))*0.1
    yv2=filX['qres'][0:tfin]
    xv2 = np.arange(0,len(yv2))*0.1

    return xv1,yv1,xv2,yv2


def updatedF7():
    tstart=10000
    noise,old,new,oerr,nerr=[],[],[],[],[]
    oldC,newC,oerrC,nerrC=[],[],[],[]
    for na in [11,12,13]:
	fileN = "analysis/random/simM_16x16_shot_n"+str(na)+"_s3.txt"
	filX=pd.read_csv(fileN)#.dropna()
	fileN = "analysis/random/simM_16x16_shot_n"+str(na)+"_s0.txt"
	fil2=pd.read_csv(fileN)#.dropna()
	noise+=[na]
	old +=[np.mean(filX['Sim'])]
	oerr +=[np.std(filX['Sim'])]
	new +=[np.mean(fil2['Sim'])]
	nerr +=[np.std(fil2['Sim'])]

	fileN = "analysis/random/Cont_16x16_shot_n"+str(na)+"_s3.txt"
	filX=pd.read_csv(fileN)#.dropna()
	fileN = "analysis/random/Cont_16x16_shot_n"+str(na)+"_s0.txt"
	fil2=pd.read_csv(fileN)#.dropna()

	oldC+=[np.mean(filX['contOpp']/5.12)]
	oerrC+=[np.std(filX['contOpp']/5.12)]
	newC+=[np.mean(fil2['contOpp']/5.12)]
	nerrC+=[np.std(fil2['contOpp']/5.12)]

    return noise, old,new,oerr,nerr,oldC,newC,oerrC,nerrC

	

###########################################################
###########################################################
########### Get data sets #################################
###########################################################
###########################################################
'''

pattI,pattF,dist,cc,time = getSpecialPatt('p0')
fileo  = open("si_data/fig2_a1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_a2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p1')
fileo  = open("si_data/fig2_b1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_b2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p2')
fileo  = open("si_data/fig2_c1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_c2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p3')
fileo  = open("si_data/fig2_d1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_d2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p4')
fileo  = open("si_data/fig2_e1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_e2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p5')
fileo  = open("si_data/fig2_f1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_f2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getSpecialPatt('p6')
fileo  = open("si_data/fig2_g1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig2_g2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
	fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()



## Figure 3
pattI,pattF,dist,cc,time = getMistakesPatt('m3','shot',0)
fileo  = open("si_data/fig3_a1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_a2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getMistakesPatt('m26','shot',0)
fileo  = open("si_data/fig3_b1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_b2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getMistakesPatt('m128','shot',0)
fileo  = open("si_data/fig3_c1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_c2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getMistakesPatt('m192','shot',0)
fileo  = open("si_data/fig3_d1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_d2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getMistakesPatt('m243','shot',0)
fileo  = open("si_data/fig3_e1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_e2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getMistakesPatt('m253','shot',0)
fileo  = open("si_data/fig3_f1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig3_f2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()


### get Figure 4
pattI,pattF,dist,cc,time = getDevPatt('d1','shot',0)
fileo  = open("si_data/fig4_a1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_a2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getDevPatt('d10','shot',0)
fileo  = open("si_data/fig4_b1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_b2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getDevPatt('d25','shot',0)
fileo  = open("si_data/fig4_c1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_c2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getDevPatt('d100','shot',0)
fileo  = open("si_data/fig4_d1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_d2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getDevPatt('d1000','white',0)
fileo  = open("si_data/fig4_e1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_e2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

pattI,pattF,dist,cc,time = getDevPatt('d2000','shot',0)
fileo  = open("si_data/fig4_f1.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattI)):
        for j in range(len(pattI[i])):
                fileo.write("%s,%s,%s,%s\n" %(pattI[i][j],pattF[i][j],i,j))
fileo.close()
fileo  = open("si_data/fig4_f2.txt",'w')
fileo.write("Sim,Cont,time\n")
for i in range(len(dist)):
        fileo.write("%s,%s,%s\n" %(dist[i],cc[i],time[i]))
fileo.close()

#### Figure 5
[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',0)
fileo  = open("si_data/fig5a0.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',7)
fileo  = open("si_data/fig5b0.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',13)
fileo  = open("si_data/fig5c0.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('shot',20)
fileo  = open("si_data/fig5d0.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',0)
fileo  = open("si_data/fig5a1.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',70)
fileo  = open("si_data/fig5b1.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',130)
fileo  = open("si_data/fig5c1.txt",'w')
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

[N,D,U,xlim,ylim,x_minima,y_minima,ntr,dtr,nts,dts]=get_pseudo('white',200)
fileo  = open("si_data/fig5d1.txt",'w')
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

### get Figure 6
xa,ya,yea,xb,yb,yeb,ya1,yb1,yea1,yeb1=getAvgs_rand()
fileo  = open("si_data/fig6.txt",'w')
fileo.write("xwhite,avgR_white,avgS_white,errR_white,errS_white,xshot,avgR_shot,avgS_shot,errR_shot,errS_shot\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],ya1[i],yea[i],yea1[i],xb[i],yb[i],yb1[i],yeb[i],yeb1[i]))
fileo.close()

### get Figure 7
xa,ya,yea,xb,yb,yeb,ya1,yb1,yea1,yeb1=getAvgs_check()
fileo  = open("si_data/fig7.txt",'w')
fileo.write("xwhite,avgR_white,avgS_white,errR_white,errS_white,xshot,avgR_shot,avgS_shot,errR_shot,errS_shot\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],ya1[i],yea[i],yea1[i],xb[i],yb[i],yb1[i],yeb[i],yeb1[i]))
fileo.close()

### get Figure 8
xa,ya,yea,xb,yb,yeb=getSimSP()
fileo  = open("si_data/fig8_check.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

xa,ya,yea,xb,yb,yeb=getSim_rand()
fileo  = open("si_data/fig8_rand.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()


#### get Figure 9
Dir='mistakes'
ext='m13'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig9_1_m13.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig9_2_m13.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='m64'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig9_1_m64.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig9_2_m64.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='m230'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig9_1_m230.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig9_2_m230.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()


#### get Figure 10
Dir='dev'
ext='d50'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig10_1_d50.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig10_2_d50.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='d1500'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig10_1_d1500.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig10_2_d1500.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='d5000'
xa,ya,yea,xb,yb,yeb=getSim(Dir,ext)
fileo  = open("si_data/fig10_1_d5000.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts(Dir,ext)
fileo  = open("si_data/fig10_2_d5000.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()


#### get Figure 11
ext='4x4'
scale=2*4*4/100.
xa,ya,yea,xb,yb,yeb=getSim_lattice(ext)
fileo  = open("si_data/fig11_1_4x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts_lattice(ext,scale)
fileo  = open("si_data/fig11_2_4x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='8x8'
scale=2*8*8/100.
xa,ya,yea,xb,yb,yeb=getSim_lattice(ext)
fileo  = open("si_data/fig11_1_8x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts_lattice(ext,scale)
fileo  = open("si_data/fig11_2_8x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

ext='32x32'
scale=2*32*32/100.
xa,ya,yea,xb,yb,yeb=getSim_lattice(ext)
fileo  = open("si_data/fig11_1_32x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()
xa,ya,yea,xb,yb,yeb=getContacts_lattice(ext,scale)
fileo  = open("si_data/fig11_2_32x.txt",'w')
fileo.write("xwhite,ywhite,yerrW,xshot,yshot,yerrS\n")
for i in range(len(xa)):
      fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()


'''



####### Figure 12

## GET FIGURE 4c
[a,b,c,std]=getET_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_white_n130_s0.txt").dropna())
fileo  = open("si_data/fig12a.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getET_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_shot_n13_s0.txt").dropna())
fileo  = open("si_data/fig12b.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getET_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_white_n200_s0.txt").dropna())
fileo  = open("si_data/fig12c.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getET_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_shot_n20_s0.txt").dropna())
fileo  = open("si_data/fig12d.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()


####### Figure 13
[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_white_n130_s0.txt").dropna())
fileo  = open("si_data/fig13a.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_shot_n13_s0.txt").dropna())
fileo  = open("si_data/fig13b.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_white_n200_s0.txt").dropna())
fileo  = open("si_data/fig13c.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effS2R_16x16_shot_n20_s0.txt").dropna())
fileo  = open("si_data/fig13d.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

####### Figure 14
[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effR2S_16x16_white_n130_s0.txt").dropna())
fileo  = open("si_data/fig14a.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effR2S_16x16_shot_n13_s0.txt").dropna())
fileo  = open("si_data/fig14b.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effR2S_16x16_white_n200_s0.txt").dropna())
fileo  = open("si_data/fig14c.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

[a,b,c,std]=getETvSR_avg(pd.read_csv("analysis/random/et_eq_effR2S_16x16_shot_n20_s0.txt").dropna())
fileo  = open("si_data/fig14d.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
        fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()


exit()


############# Figure 15
x,n0w,n0s,n50,n130,n200,n5,n13,n20=getSim_checkFull_noise()
fileo  = open("si_data/fig15_checkS.txt",'w')
fileo.write("x,n0w,n0s,n50,n130,n200,n5,n13,n20\n")
for i in range(len(x)):
      fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(x[i],n0w[i],n0s[i],n50[i],n130[i],n200[i],n5[i],n13[i],n20[i]))
fileo.close()
x,n0w,n0s,n50,n130,n200,n5,n13,n20=getCont_checkFull_noise()
fileo  = open("si_data/fig15_checkC.txt",'w')
fileo.write("x,n0w,n0s,n50,n130,n200,n5,n13,n20\n")
for i in range(len(x)):
      fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(x[i],n0w[i],n0s[i],n50[i],n130[i],n200[i],n5[i],n13[i],n20[i]))
fileo.close()


#### Figure 16
xv1,yv1,xv2,yv2=qcorr_rand('white',50)
fileo  = open("si_data/fig16_a50.txt",'w')
fileo.write("time,q\n")
for i in range(len(xv2)):
      fileo.write("%s,%s\n" %(xv2[i],yv2[i]))
fileo.close()

xv1,yv1,xv2,yv2=qcorr_rand('white',130)
fileo  = open("si_data/fig16_a130.txt",'w')
fileo.write("time,q\n")
for i in range(len(xv2)):
      fileo.write("%s,%s\n" %(xv2[i],yv2[i]))
fileo.close()

xv1,yv1,xv2,yv2=qcorr_rand('shot',5)
fileo  = open("si_data/fig16_b5.txt",'w')
fileo.write("time,q\n")
for i in range(len(xv2)):
      fileo.write("%s,%s\n" %(xv2[i],yv2[i]))
fileo.close()

xv1,yv1,xv2,yv2=qcorr_rand('shot',13)
fileo  = open("si_data/fig16_b13.txt",'w')
fileo.write("time,q\n")
for i in range(len(xv2)):
      fileo.write("%s,%s\n" %(xv2[i],yv2[i]))
fileo.close()

#### Updated Fig S7


noise, old,new,oerr,nerr,oldC,newC,oerrC,nerrC=updatedF7()
fileo  = open("si_data/figu7.txt",'w')
fileo.write("noise,simOld,simnew,simOerr,simNerr,contOld,contnew,contOerr,contNerr\n")
for i in range(len(noise)):
      fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(noise[i],old[i],new[i],oerr[i],nerr[i],oldC[i],newC[i],oerrC[i],nerrC[i]))
fileo.close()

[dN,dD,dI] = get_data("data/random/traj_16x16x1_white_n200.dat",tstart=0)
file1 = open("si_data/fig3_a_w2.txt",'w')
for i in range(len(dD[0])):
	file1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(dD[0][i][0][0],dD[0][i][0][1],dD[0][i][0][2],dD[0][i][0][3],dD[0][i][0][4],dD[0][i][0][5],dD[0][i][0][6],dD[0][i][0][7],dD[0][i][0][8],dD[0][i][0][9],dD[0][i][0][10],dD[0][i][0][11],dD[0][i][0][12],dD[0][i][0][13],dD[0][i][0][14],dD[0][i][0][15]))
file1.close()


[dN,dD,dI] = get_data("data/random/traj_16x16x1_white_n0.dat",tstart=0)
file1 = open("si_data/fig3_a_w0.txt",'w')
for i in range(len(dD[0])):
	file1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(dD[0][i][0][0],dD[0][i][0][1],dD[0][i][0][2],dD[0][i][0][3],dD[0][i][0][4],dD[0][i][0][5],dD[0][i][0][6],dD[0][i][0][7],dD[0][i][0][8],dD[0][i][0][9],dD[0][i][0][10],dD[0][i][0][11],dD[0][i][0][12],dD[0][i][0][13],dD[0][i][0][14],dD[0][i][0][15]))

file1.close()

[dN,dD,dI] = get_data("data/random/traj_16x16x1_shot_n20.dat",tstart=0)
file1 = open("si_data/fig3_a_s2.txt",'w')
for i in range(len(dD[0])):
	file1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(dD[0][i][0][0],dD[0][i][0][1],dD[0][i][0][2],dD[0][i][0][3],dD[0][i][0][4],dD[0][i][0][5],dD[0][i][0][6],dD[0][i][0][7],dD[0][i][0][8],dD[0][i][0][9],dD[0][i][0][10],dD[0][i][0][11],dD[0][i][0][12],dD[0][i][0][13],dD[0][i][0][14],dD[0][i][0][15]))
file1.close()

[dN,dD,dI] = get_data("data/random/traj_16x16x1_shot_n0.dat",tstart=0)
file1 = open("si_data/fig3_a_s0.txt",'w')
for i in range(len(dD[0])):
	file1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(dD[0][i][0][0],dD[0][i][0][1],dD[0][i][0][2],dD[0][i][0][3],dD[0][i][0][4],dD[0][i][0][5],dD[0][i][0][6],dD[0][i][0][7],dD[0][i][0][8],dD[0][i][0][9],dD[0][i][0][10],dD[0][i][0][11],dD[0][i][0][12],dD[0][i][0][13],dD[0][i][0][14],dD[0][i][0][15]))
file1.close()



for sim in range(20):
	file1 = open("si_data/figS_ind_s"+str(sim)+".txt",'w')
	file1.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
	data= getSim_rand_ind(sim)
	for i in range(len(data[0])):
		file1.write("%s,%s,%s,%s,%s,%s\n" %(data[0][i],data[1][i],data[2][i],data[3][i],data[4][i],data[5][i]))
	file1.close()

	file1 = open("si_data/figC_ind_s"+str(sim)+".txt",'w')
	file1.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
	data= getCont_rand_ind(sim)
	for i in range(len(data[0])):
		file1.write("%s,%s,%s,%s,%s,%s\n" %(data[0][i],data[1][i],data[2][i],data[3][i],data[4][i],data[5][i]))
	file1.close()


