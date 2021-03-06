## Written by Madeline Galbraith
## Last edited: July 2022

import numpy as np
import pandas as pd
from auxFunctions import get_data
from auxFunctions import get_thresholds
from plot_script_functions import get_dataVersion2
from plot_script_functions import pseudo_potential
from plot_script_functions import find_minima

#############################
#############################
#####  Figure 2 data ########
#############################
#############################
def getCont_mist_det():
    tstart=10000
    xvals,yvals,yerr=[],[],[]
    for extV in [3,13,26,64,128,192,230,243,253]:
    	tmp=[]
        ext = 'm'+str(extV)
        for sim in range(20):
            try:
                fileN = "data/mistakes/Cont_16x16_"+str(ext)+"_"+str(sim)+".txt"
                filX=pd.read_csv(fileN).dropna()
                tmp+=list(filX['contOpp'].values[tstart:])
            except:
                count=0

	tmp = np.array(tmp)/5.12
        xvals+=[float(ext.replace('m',''))/256.*100.]
        yvals+=[np.mean(tmp)]
        yerr+=[np.std(tmp)]
    return  xvals,yvals,yerr

def getCont_dev_det():
    tstart=10000
    xvals,yvals,yerr=[],[],[]
    for extV in [0,5,10,25,100,1000,1500,2000,5000]:
    	tmp=[]
        #for extV in [0,1,2,3,4,5,10,15,25,50,75,100,200,500,1000,1500,2000,5000]:
        ext = 'd'+str(extV)
        for sim in range(20):
            try:
                fileN = "data/dev/Cont_16x16_"+str(ext)+"_"+str(sim)+".txt"
                filX=pd.read_csv(fileN)#.dropna()
                tmp=tmp+list(filX['contOpp'].values[tstart:])
            except:
                count=0

	print np.mean(tmp),ext
	tmp = np.array(tmp)/5.12
	print np.mean(tmp),ext
        xvals+=[int(ext.replace("d",""))/1000.]
        yvals+=[np.mean(tmp)]
        yerr+=[np.std(tmp)]
    return  xvals,yvals,yerr

def getSim_mist_det():
    tstart=10000
    xvals,yvals,yerr=[],[],[]
    for extV in [3,13,26,64,128,192,230,243,253]:
    	tmp=[]
        ext = 'm'+str(extV)
        for sim in range(20):
            try:
                fileN = "data/mistakes/simM_16x16_"+str(ext)+"_"+str(sim)+".txt"
                filX=pd.read_csv(fileN).dropna()
                tmp+=list(filX['Sim'].values[tstart:])
            except:
                count=0
        xvals+=[float(ext.replace('m',''))/256.*100.]
        yvals+=[np.mean(tmp)]
        yerr+=[np.std(tmp)]
    return  xvals,yvals,yerr

def getSim_dev_det():
    tstart=10000
    xvals,yvals,yerr=[],[],[]
    for extV in [0,5,10,25,100,1000,1500,2000,5000]:
    	tmp=[]
        #for extV in [0,1,2,3,4,5,10,15,25,50,75,100,200,500,1000,1500,2000,5000]:
        ext = 'd'+str(extV)
        for sim in range(20):
            try:
                fileN = "data/dev/simM_16x16_"+str(ext)+"_"+str(sim)+".txt"
                filX=pd.read_csv(fileN).dropna()
                tmp=tmp+list(filX['Sim'].values[tstart:])
            except:
                count=0

        xvals+=[int(ext.replace("d",""))/1000.]
        yvals+=[np.mean(tmp)]
        yerr+=[np.std(tmp)]
    return  xvals,yvals,yerr

def getSim_randFull_det():
    fileN = "data/random/simM_16x16_shot_n0_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    ind= np.argwhere(np.isnan(filX['Sim']))[:,0]
    filX['Sim'][ind]=0
    x =np.arange(0,len(filX['Sim']))*0.1
    return x,filX['Sim']
def getCont_randFull_det():
    fileN = "data/random/Cont_16x16_shot_n0_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    ind= np.argwhere(np.isnan(filX['contOpp']))[:,0]
    filX['contOpp'][ind]=0
    x =np.arange(0,len(filX['contOpp']))*0.1
    return x,filX['contOpp']


def getSim_checkFull_det():
    #fileN = "data/patt/simM_16x16_shot_n0_p8_s0.txt"
    fileN = "data/patt/simM_16x16_shot_n0_p7_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    ind= np.argwhere(np.isnan(filX['Sim']))[:,0]
    filX['Sim'][ind]=0
    x =np.arange(0,len(filX['Sim']))*0.1
    return x,filX['Sim']

def getCont_checkFull_det():
    #fileN = "data/patt/Cont_16x16_shot_n0_p8_s0.txt"
    fileN = "data/patt/Cont_16x16_shot_n0_p7_s0.txt"
    filX=pd.read_csv(fileN)#.dropna()
    ind= np.argwhere(np.isnan(filX['contOpp']))[:,0]
    filX['contOpp'][ind]=0
    x =np.arange(0,len(filX['contOpp']))*0.1
    return x,filX['contOpp']


def getCheckPatt():
    #[dN,dD,dI] = get_data("data/patt/traj_16x16x1_shot_n0_p8_s0.dat",tstart=0)
    [dN,dD,dI] = get_data("data/patt/traj_16x16x1_shot_n0_p7_s0.dat",tstart=0)
    kk = dN.keys()[0]
    
    if False:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)

        iD = dtr+ (dts-dtr)/2.
        iN = nts+ (ntr-nts)/2.
        pattI = (dD[kk][0]>iD)*(dN[kk][0]<iN)*1.+(dD[kk][0]<iD)*(dN[kk][0]>iN)*2.
        pattF = (dD[kk][-1]>iD)*(dN[kk][-1]<iN)*1.+(dD[kk][-1]<iD)*(dN[kk][-1]>iN)*2.
    if True:
        pattI = dD[kk][0]
        pattF = dD[kk][-1]        
    return pattI,pattF

def getrandPatt():
    [dN,dD,dI] = get_data("data/random/traj_16x16x1_white_n0.dat",tstart=0)
    kk = dN.keys()[0]
    if False:
        [nr,dr,ns,ds,ntr,dtr,nts,dts] = get_thresholds(dN,dD)

        iD = dtr+ (dts-dtr)/2.
        iN = nts+ (ntr-nts)/2.
        pattI = (dD[kk][0]>iD)*(dN[kk][0]<iN)*1.+(dD[kk][0]<iD)*(dN[kk][0]>iN)*2.
        pattF = (dD[kk][-1]>iD)*(dN[kk][-1]<iN)*1.+(dD[kk][-1]<iD)*(dN[kk][-1]>iN)*2.
    if True:
        pattI = dD[kk][0]
        pattF = dD[kk][-1] 
    return pattI,pattF

'''
## GET FIGURE 2A
pattIr,pattFr = getrandPatt()
fileo  = open("main_data/fig2a.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattIr)):
	for j in range(len(pattIr[i])):
		fileo.write("%s,%s,%s,%s\n" %(pattIr[i][j],pattFr[i][j],i,j))
fileo.close()

## GET FIGURE 2B
pattIc,pattFc = getCheckPatt()
fileo  = open("main_data/fig2b.txt",'w')
fileo.write("Initial,Final,row,column\n")
for i in range(len(pattIc)):
	for j in range(len(pattIc[i])):
		fileo.write("%s,%s,%s,%s\n" %(pattIc[i][j],pattFc[i][j],i,j))
fileo.close()

## GET FIGURE 2C
xc,yc = getSim_checkFull_det()
xr,yr = getSim_randFull_det()
fileo  = open("main_data/fig2c.txt",'w')
fileo.write("xcheck,ycheck,xrand,yrand\n")
for i in range(len(xc)):
	fileo.write("%s,%s,%s,%s\n" %(xc[i],yc[i],xr[i],yr[i]))
fileo.close()

## GET FIGURE 2C--cont
xc,yc = getCont_checkFull_det()
xr,yr = getCont_randFull_det()
fileo  = open("main_data/fig2c_cont.txt",'w')
fileo.write("xcheck,ycheck,xrand,yrand\n")
for i in range(len(xc)):
	fileo.write("%s,%s,%s,%s\n" %(xc[i],yc[i],xr[i],yr[i]))
fileo.close()

## GET FIGURE 2D
x,y,yerr = getSim_mist_det()
fileo  = open("main_data/fig2d-opt2.txt",'w')
fileo.write("x,y,yerr\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s\n" %(x[i],y[i],yerr[i]))
fileo.close()

## GET FIGURE 2E
x,y,yerr = getSim_dev_det()
fileo  = open("main_data/fig2e-opt2.txt",'w')
fileo.write("x,y,yerr\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s\n" %(x[i],y[i],yerr[i]))
fileo.close()


## GET FIGURE 2D
x,y,yerr = getCont_mist_det()
fileo  = open("main_data/fig2d.txt",'w')
fileo.write("x,y,yerr\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s\n" %(x[i],y[i],yerr[i]))
fileo.close()

## GET FIGURE 2E
x,y,yerr = getCont_dev_det()
fileo  = open("main_data/fig2e.txt",'w')
fileo.write("x,y,yerr\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s\n" %(x[i],y[i],yerr[i]))
fileo.close()
'''

#############################
#############################
#####  Figure 3 data ########
#############################
#############################

def getAvgs_rand():
    tstart=10000
    tfin=100000
    xa,xb=[],[]
    ya1,ya2,yb1,yb2=[],[],[],[]
    yera1,yera2,yerb1,yerb2=[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp,tmp1=[],[]
            sim=0
	    for sim in range(20):
                try:
                        fileN = "data/random/avgStates_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                        filX=pd.read_csv(fileN).dropna()
                        #tmp+=list(filX['avgR'].values[tstart:tfin]+filX['avglR'].values[tstart:tfin])
                        tmp+=list(filX['avgER'].values[tstart:tfin])
                        #tmp1+=list(filX['avgS'].values[tstart:tfin]+filX['avglS'].values[tstart:tfin])
                        tmp1+=list(filX['avgES'].values[tstart:tfin])
                except:
                    count=0
            tmp = np.array(tmp)
            tmp1 = np.array(tmp1)
            if nt=='white':
                xa+=[na]
                ya1+=[np.mean(tmp*100.)]
                ya2+=[np.mean(tmp1*100.)]
                yera1+=[np.std(tmp*100.)]
                yera2+=[np.std(tmp1*100.)]
            elif nt=='shot':
                xb+=[na]
                yb1+=[np.mean(tmp*100.)]
                yb2+=[np.mean(tmp1*100.)]
                yerb1+=[np.std(tmp*100.)]
                yerb2+=[np.std(tmp1*100.)]
    return xa,ya1,yera1,xb,yb1,yerb1,ya2,yb2,yera2,yerb2

def getContacts_rand():
    tstart=10000
    tfin=100000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "data/random/Cont_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['contOpp'].values[tstart:tfin])
                except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp/5.12)]
                yera+=[np.std(tmp/5.12)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp/5.12)]
                yerb+=[np.std(tmp/5.12)]
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
                    fileN = "data/random/simM_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
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
                    fileN = "data/patt/simM_16x16_"+nt+"_n"+str(na)+"_p8_s"+str(sim)+".txt"
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

def getContactsSP():
    tstart=10000
    xa,xb,ya,yb,yera,yerb=[],[],[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp=[]
            for sim in range(20):
                try:
                    fileN = "data/patt/Cont_16x16_"+nt+"_n"+str(na)+"_p8_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN).dropna()
                    tmp+=list(filX['contOpp'].values[tstart:])
                except:
                    count=0
            tmp = np.array(tmp)
            if nt=='white':
                xa+=[na]
                ya+=[np.mean(tmp/5.12)]
                yera+=[np.std(tmp/5.12)]
            elif nt=='shot':
                xb+=[na]
                yb+=[np.mean(tmp/5.12)]
                yerb+=[np.std(tmp/5.12)]
    return xa,ya,yera,xb,yb,yerb

def getAvgsSP():
    tstart=10000
    xa,xb=[],[]    
    ya1,ya2,yb1,yb2=[],[],[],[]
    yera1,yera2,yerb1,yerb2=[],[],[],[]
    nAmp={'shot':np.arange(0,21,1),'white':np.arange(0,210,10)}
    for nt in ['shot','white']:
        for na in nAmp[nt]:
            tmp,tmp1,tmp2,tmp3=[],[],[],[]
            for sim in range(20):
                try:
                    fileN = "data/patt/avgStates_16x16_"+nt+"_n"+str(na)+"_p8_s"+str(sim)+".txt"
                    filX=pd.read_csv(fileN).dropna()
                    #tmp+=list(filX['avgR'].values[tstart:]+filX['avglR'].values[tstart:])
                    tmp+=list(filX['avgER'].values[tstart:])
                    #tmp1+=list(filX['avgS'].values[tstart:]+filX['avglS'].values[tstart:])
                    tmp1+=list(filX['avgES'].values[tstart:])
                except:
                    count=0
            tmp = np.array(tmp)
            tmp1 = np.array(tmp1)
            if nt=='white':
                xa+=[na]
                ya1+=[np.mean(tmp*100.)]
                ya2+=[np.mean(tmp1*100.)]
                yera1+=[np.std(tmp*100.)]
                yera2+=[np.std(tmp1*100.)]
            elif nt=='shot':
                xb+=[na]
                yb1+=[np.mean(tmp*100.)]
                yb2+=[np.mean(tmp1*100.)]
                yerb1+=[np.std(tmp*100.)]
                yerb2+=[np.std(tmp1*100.)]
    return xa,ya1,yera1,xb,yb1,yerb1,ya2,yb2,yera2,yerb2

def get_pseudo(nt,na):
    # parameters for plotting
    sim = 20
    #sim_points = 10000
    sim_points = 20000#150000
    #n_eq = 1000
    n_eq = 10000
    bins = np.logspace(np.log10(0.01), np.log10(100000), num=61, base=10)
    ncell = 16

    N, D = get_dataVersion2('data_ND_mcT_trajs_Oct2020/traj_16x16x1_'+str(nt)+'_n'+str(na)+'.dat',10000,25000,0)
    U = pseudo_potential(N, D, bins)
    
    x = np.zeros(bins.size - 1)
    for i in range(x.size):
        x[i] = (bins[i] + bins[i + 1]) / 2.
        
    notch=bins[0:-1]
    delta=bins[0:-1]
    xint=[10, 30000]
    yint=[0.01, 10000]

    # find minima
    i_s, j_s, i_r, j_r = find_minima(U, bins, x, x)
    
    return [notch,delta,np.transpose(U),xint,yint,[bins[i_s], bins[i_r]], [bins[j_s], bins[j_r]]]


## GET FIGURE 3A
'''
[N,D,U,xlim,ylim,x_minima,y_minima]=get_pseudo('shot',0)

fileo  = open("main_data/fig3a.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima\n")
count=0
for i in range(len(N)):
	for j in range(len(D)):
		if count==0:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0]))
		elif count==1:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1]))
		else:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1))
		count+=1
fileo.close()

## GET FIGURE 3B
[N,D,U,xlim,ylim,x_minima,y_minima]=get_pseudo('shot',10)
fileo  = open("main_data/fig3b.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima\n")
count=0
for i in range(len(N)):
	for j in range(len(D)):
		if count==0:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0]))
		elif count==1:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1]))
		else:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1))
		count+=1
fileo.close()

## GET FIGURE 3C
[N,D,U,xlim,ylim,x_minima,y_minima]=get_pseudo('shot',20)
fileo  = open("main_data/fig3c.txt",'w')
fileo.write("N,D,U,xlim,ylim,xminima,yminima\n")
count=0
for i in range(len(N)):
	for j in range(len(D)):
		if count==0:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[0],ylim[0],x_minima[0],y_minima[0]))
		elif count==1:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],xlim[1],ylim[1],x_minima[1],y_minima[1]))
		else:
			fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(N[i],D[j],U[i][j],-1,-1,-1,-1))
		count+=1
fileo.close()

## GET FIGURE 3D
xa,ya,yea,xb,yb,yeb=getContactsSP()
fileo  = open("main_data/fig3d.txt",'w')
fileo.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
for i in range(len(xa)):
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

## GET FIGURE 3E
xa,ya,yea,xb,yb,yeb=getContacts_rand()
fileo  = open("main_data/fig3e.txt",'w')
fileo.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
for i in range(len(xa)):
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()


## GET FIGURE 3D
xa,ya,yea,xb,yb,yeb=getSimSP()
fileo  = open("main_data/fig3d2.txt",'w')
fileo.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
for i in range(len(xa)):
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

## GET FIGURE 3E
xa,ya,yea,xb,yb,yeb=getSim_rand()
fileo  = open("main_data/fig3e2.txt",'w')
fileo.write("xwhite,ywhite,yewhite,xshot,yshot,yeshot\n")
for i in range(len(xa)):
	fileo.write("%s,%s,%s,%s,%s,%s\n" %(xa[i],ya[i],yea[i],xb[i],yb[i],yeb[i]))
fileo.close()

'''

#############################
#############################
#####  Figure 4 data ########
#############################
#############################

def getCont_randFull_noise():
    tfin=40000
    fileN = "data/random/Cont_16x16_white_n0_s0.txt"
    filX=pd.read_csv(fileN)
    n0w=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_white_n50_s0.txt"
    filX=pd.read_csv(fileN)
    n50=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_white_n130_s0.txt"
    filX=pd.read_csv(fileN)
    n130=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_white_n200_s0.txt"
    filX=pd.read_csv(fileN)
    n200=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_white_n0_s0.txt"
    filX=pd.read_csv(fileN)
    n0s=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_shot_n5_s0.txt"
    filX=pd.read_csv(fileN)
    n5=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_shot_n13_s0.txt"
    filX=pd.read_csv(fileN)
    n13=filX['contOpp'][:tfin]/5.12
    fileN = "data/random/Cont_16x16_shot_n20_s0.txt"
    filX=pd.read_csv(fileN)
    n20=filX['contOpp'][:tfin]/5.12

    x =np.arange(0,len(filX['contOpp'][:tfin]))*0.1
    return x,n0w,n0s,n50,n130,n200,n5,n13,n20

def getSim_randFull_noise():
    tfin=40000
    fileN = "data/random/simM_16x16_white_n0_s0.txt"
    filX=pd.read_csv(fileN)
    n0w=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_white_n50_s0.txt"
    filX=pd.read_csv(fileN)
    n50=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_white_n130_s0.txt"
    filX=pd.read_csv(fileN)
    n130=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_white_n200_s0.txt"
    filX=pd.read_csv(fileN)
    n200=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_white_n0_s0.txt"
    filX=pd.read_csv(fileN)
    n0s=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_shot_n5_s0.txt"
    filX=pd.read_csv(fileN)
    n5=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_shot_n13_s0.txt"
    filX=pd.read_csv(fileN)
    n13=filX['Sim'][:tfin]
    fileN = "data/random/simM_16x16_shot_n20_s0.txt"
    filX=pd.read_csv(fileN)
    n20=filX['Sim'][:tfin]

    x =np.arange(0,len(filX['Sim'][:tfin]))*0.1
    return x,n0w,n0s,n50,n130,n200,n5,n13,n20

def getET_avg(df):
    x,y = np.round(df['D'],0),np.round(df['N'],-1)
    yval=(np.unique(y))
    xval=(np.unique(x))

    a,b,c=[],[],[]
    std=[]

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

## GET FIGURE 4a
[a,b,c,std]=getET_avg( pd.read_csv("data/random/et_eq_effR2S_16x16_white_n130_s0.txt").dropna())
fileo  = open("main_data/fig4a.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
	fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

## GET FIGURE 4b
[a,b,c,std]=getET_avg( pd.read_csv("data/random/et_eq_effR2S_16x16_shot_n13_s0.txt").dropna())
fileo  = open("main_data/fig4b.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
	fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

## GET FIGURE 4c
[a,b,c,std]=getET_avg( pd.read_csv("data/random/et_eq_effR2S_16x16_white_n200_s0.txt").dropna())
fileo  = open("main_data/fig4c.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
	fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

## GET FIGURE 4d
[a,b,c,std]=getET_avg( pd.read_csv("data/random/et_eq_effR2S_16x16_shot_n20_s0.txt").dropna())
fileo  = open("main_data/fig4d.txt",'w')
fileo.write("x,y,color\n")
for i in range(len(a)):
	fileo.write("%s,%s,%s\n" %(a[i],b[i],c[i]))
fileo.close()

## GET FIGURE 4ef
[x,n0w,n0s,n50,n130,n200,n5,n13,n20] = getSim_randFull_noise()
fileo  = open("main_data/fig4ef-opt2.txt",'w')
fileo.write("time,white0,white50,white130,white200,shot0,shot5,shot13,shot20\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(x[i],n0w[i],n50[i],n130[i],n200[i],n0s[i],n5[i],n13[i],n20[i]))
fileo.close()

## GET FIGURE 4ef
[x,n0w,n0s,n50,n130,n200,n5,n13,n20] = getCont_randFull_noise()
fileo  = open("main_data/fig4ef.txt",'w')
fileo.write("time,white0,white50,white130,white200,shot0,shot5,shot13,shot20\n")
for i in range(len(x)):
	fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(x[i],n0w[i],n50[i],n130[i],n200[i],n0s[i],n5[i],n13[i],n20[i]))
fileo.close()

