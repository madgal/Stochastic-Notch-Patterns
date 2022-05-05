import numpy as np
import pandas as pd
from scipy.stats import pearsonr


### The functions include ####
'''
def convert2timeSeries(dat,size,tstart):
def get_thresholds(dN,dD):
def readFile(fn,dirct,tstart):
def find_thresholds(x, U, i_s, j_s, i_r, j_r):
def find_minima(dist, xbins, notch, delta):
def pseudo_potential(N, D, bins):
def get_data(filen,dirct,tstart=1000):
def getCorr(X):
def getSimStates(dD,dN,tD,tN):
def getET_states(states,dN,dD,tfin=-1):
def get_results1(ddN,ddD,thresholds):
def get_results2(ddN,ddD,thresholds):
def get_results3(ddN,ddD,thresholds):
def get_results7(ddN,ddD,thresholds):
def getCont(states):
'''
	
__thresholds__={}

def get_thresholds(dN,dD):
    ## based off code written by Federico for finding the minima and getting the threshold at 0.1 of U

    key = list(dN.keys())[0]
    ## select the correct data from the file
    N= dN[key]
    D= dD[key]

	
    ## create the bins for the pseudopotential
    bins = np.logspace(np.log10(0.1),np.log10(100000),num=121,base=10)
    x = np.zeros(bins.size-1)
    for i in range(x.size):
        x[i] = (bins[i]+bins[i+1])/2.

    ## generate U and get the key information
    U = pseudo_potential(N,D,bins)
    i_s,j_s,i_r,j_r = find_minima(U,bins,x,x)
    ntr,dtr,nts,dts = find_thresholds(x,U,i_s,j_s,i_r,j_r)
    nr = bins[i_s]## Notch threshold for Receiver
    dr = bins[j_s]## Delta threshold for Receiver
    ns = bins[i_r]## Notch threshold for Sender
    ds = bins[j_r]## Delta threshold for Sender

    return [nr,dr,ns,ds,ntr,dtr,nts,dts]

def convert2timeSeries(dat,size,tstart):
	#### functions for reading and analyzing data
        ## takes the column from the trajectory files and
        ## converts it into a time series with each step as an NxN lattice

        dat=np.array(dat)
        for row in range(len(dat)):
                dat[row] = dat[row].replace("[","")
                dat[row] = dat[row].replace("]","")
                dat[row] = np.array(dat[row].split(","))
                dat[row] = dat[row].astype(float)
                dat[row] = np.reshape(dat[row],(size,size))
        dat = dat[tstart:]
	
        return dat

def readFile(fn,tstart):
	## read the file and put it into the appropriate format
	## catch and fix any exceptions that may occur during reading 

        df=  pd.read_csv(filepath_or_buffer=fn,header=0).values

        filex=(fn.split(".")[0]).split("_")
        lattice_size = filex[1]
        lattice_size = int(lattice_size.split("x")[0])

        ## data for N,D, and I respectively in the form
        ## sn[simulationNumber] = [ [NxN lattice at t=0],[NxN lattice at t=1]...,[NxN lattice at t_final]]
        sn,sd,si={},{},{}
        Nsim =np.unique(df[:,3]) ##get the number of simulations present (should be 20), can also select just a few runs

        for simv in Nsim:
		fail=False
		try:
                        ind = np.argwhere(df[:,3]==simv)[:,0] ## indices of all data for that simulation, will be consecutive
                        simv = int(simv) ## make sure it's an integer
			if len(ind)==0:
				fail=True
                        ## convert the columns indivdually
                        sn[simv] = convert2timeSeries(df[ind,2],lattice_size,tstart)
                        sd[simv] = convert2timeSeries(df[ind,0],lattice_size,tstart)
                        si[simv] = convert2timeSeries(df[ind,1],lattice_size,tstart)
		except:
			print simv," not converted to int"
			fail=True
		if fail:
			try:
                                simv = int(simv) ## make sure it's an integer
                                ind = np.argwhere(df[:,3]==simv)[:,0] ## indices of all data for that simulation, will be consecutive
                                ## convert the columns indivdually
                                sn[simv] = convert2timeSeries(df[ind,2],lattice_size,tstart)
                                sd[simv] = convert2timeSeries(df[ind,0],lattice_size,tstart)
                                si[simv] = convert2timeSeries(df[ind,1],lattice_size,tstart)
			except:
				print simv," not converted to int"
        return [sn,sd,si]

def get_data(filen,tstart=1000):
    ## return the data as a tuple of Notch, Delta, and NICD copy numbers
    [dataN,dataD,dataI] = readFile(filen,tstart)
    myData=[dataN,dataD,dataI]
    return myData

def getCorr(X):
	## get the correlation of the lattice
	## use the pearson correlation across the rows and columns, then average
        [nr,nc] = X.shape

	## pearson correlation of columns
        totalC,count=0,0
        for i in range(nc):
                corr,_ = pearsonr(X[i],X[(i+1)%nc])

                corr = (1.-corr)/2 # convert from (-1,1) to (0,1)
                totalC+=corr
                count+=1.
        totalC=totalC/count

	## pearson correlation of rows
        totalR,count=0,0
        for i in range(nr):
                corr,_ = pearsonr(X[:,i],X[:,(i+1)%nr])

                corr = (1.-corr)/2 # convert from (-1,1) to (0,1)
                totalR+=corr
                count+=1.
        totalR=totalR/count

        return (totalR+totalC)/2.

def getSimStates(dD,dN,tD,tN):
    ## dD is Delta values
    ## dN is Notch values
    ## tD is threshold values for delta
    ## tN is threshold values for notch
    ## go through all the cells in the lattice and find how they behave

    states=[]## define S or R based on whether Notch or Delta greater, all cells either S=1 or R=2
    states2=[]## define S or R based on midpoint, all cells either S=1 or R=2 
    states3=[]## define S or R based on midpoint, all cells either S=1 or R=-1
    ## find midpoint between Sender and Receiver to ensure all cells are classified as either S or R
    iD = tD['R']+ (tD['S']-tD['R'])/2. 
    iN = tN['S']+ (tN['R']-tN['S'])/2. 
    for k in range(len(dD)):
        if k==0:
            ### set the initial as either sender or receiver by breaking the thresholds in halfo
	    tmp = (dD[k]>=dN[k])*1.+(dD[k]<dN[k])*2.

        states+=[tmp]

        tmpB = (dD[k]>iD)*(dN[k]<iN)*1. +(dD[k]<iD)*(dN[k]>iN)*2.
        states2+=[tmpB]

        tmpC = (dD[k]>iD)*(dN[k]<iN)*1. -(dD[k]<iD)*(dN[k]>iN)*1.
        states3+=[tmpC]

    return [states,states2,states3]

def updateDictionary(dict,Nv,Dv,Rv,Sv):
        dict['N']+=[np.mean(tmpNs)]
        dict['D']+=[np.mean(tmpDs)]
        dict['R']+=[np.mean(numRs)]
        dict['S']+=[np.mean(numSs)]
	return dict
def updateValue(data,rv,cv,Nrow,Ncol):
    return (data[(rv-1)%Nrow][cv]+data[(rv+1)%Nrow][cv]+data[rv][(cv-1)%Ncol]+data[rv][(cv+1)%Ncol])/4.
def updateState(data,comp,Nrow,Ncol):
    return ((data[(rv-1)%Nrow][cv]==comp)*1.+(data[(rv+1)%Nrow][cv]==comp)*1.+(data[rv][(cv-1)%Ncol]==comp)*1.+(data[rv][(cv+1)%Ncol]==comp)*1.)

def getET_states(states,dN,dD,tfin=-1):
    ## first determine how long it is in S or R state
    ## save that in resS or resR
    ## then determine how long in transition state 
    ## the transition from R to S means: start adding time once it enters the R and until it enters S
    ## S to R is vice versa
    times ={}
    neighS,neighV={},{}
    for sim in states:
        times[sim]={}
        neighS[sim]={}
        neighV[sim]={}
        Nrow,Ncol= states[sim][0].shape
        for rv in range(Nrow):
            times[sim][rv]={}
            neighV[sim][rv]={}
            neighS[sim][rv]={}
            for cv in range(Ncol):
                times[sim][rv][cv]={'effS2R':[],'effR2S':[],'resS':[],'resR':[],'resIr':[],'resIs':[]}
                neighV[sim][rv][cv]={'effS2R':{'N':[],'D':[]},'effR2S':{'N':[],'D':[]},'resS':{'N':[],'D':[]},'resR':{'N':[],'D':[]},'resIr':{'N':[],'D':[]},'resIs':{'N':[],'D':[]}}
                neighS[sim][rv][cv]={'effS2R':{'S':[],'R':[]},'effR2S':{'S':[],'R':[]},'resS':{'S':[],'R':[]},'resR':{'S':[],'R':[]},'resIr':{'S':[],'R':[]},'resIs':{'S':[],'R':[]}}
                
                ### use last half of simulation unless it is very short (<20000 steps)
		if len(states[sim])>20000:
                	k=10000            
		else:
                	k=1  ## k=0 will be previous state

		## initialize previous state,prevS1
                prevS1 = states[sim][k-1][rv][cv] ## should be  1-S or 2-R
		## stTime = residency time of S or R
		## intTime = time in intermediate state (not in S or R)
		## eTime = time starting when entering R (or S) state and until it crosses to S (or R) state
                stTime,intTime,eTime=1,1,1

		### we will take the average of the neighbors values over the time of transition
		## average N,D,S,R when in S or R
                tmpNs,tmpDs,numSs,numRs=[],[],[],[]
		## average N,D,S,R when in intermediate
                tmpNi,tmpDi,numSi,numRi=[],[],[],[]
		## average N,D,S,R for entire transition
                numSe,numRe,tmpDe,tmpNe=[],[],[],[]

		## confirm the final timestep is valid or set it to be at end of simulation
		if tfin!=-1 and len(states[sim])>=tfin:
			final_step =tfin
		else:
			final_step =len(states[sim])

		## calculate times
                while (k <final_step):
                    current = states[sim][k][rv][cv]
                    if current!=prevS1:
                        if prevS1==1:# prev was S
                            times[sim][rv][cv]['resS']+=[stTime/10.]
                            neighV[sim][rv][cv]['resS'] = updateDictionary(neighV[sim][rv][cv]['resS'],tmpNs,tmpDs,numRs,numSs)
                            if current==2:## currently R and previous was S == cell transitioned
                                times[sim][rv][cv]['effS2R']+=[eTime/10.]
                                neighV[sim][rv][cv]['effS2R'] = updateDictionary(neighV[sim][rv][cv]['effS2R'],tmpNe,tmpDe,numRe,numSe)
				## reset transition variables
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            intTime=1
                            numSi,numRi,tmpDi,tmpNi=[],[],[],[]
                        elif prevS1==2:# prev was R
                            times[sim][rv][cv]['resR']+=[stTime/10.]
                            neighV[sim][rv][cv]['resR'] = updateDictionary(neighV[sim][rv][cv]['resR'],tmpNs,tmpDs,numRs,numSs)
                            if current==1:
                                times[sim][rv][cv]['effR2S']+=[eTime/10.]
                                neighV[sim][rv][cv]['effR2S'] = updateDictionary(neighV[sim][rv][cv]['effR2S'],tmpNe,tmpDe,numRe,numSe)
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            intTime=1
                            numSi,numRi,tmpDi,tmpNi=[],[],[],[]
                        elif prevS1==3:# prev was transS
                            times[sim][rv][cv]['resIs']+=[intTime/10.]
                            neighV[sim][rv][cv]['resIs'] = updateDictionary(neighV[sim][rv][cv]['resIs'],tmpNi,tmpDi,numRi,numSi)
                            if current==2:# s2r
                                times[sim][rv][cv]['effS2R']+=[eTime/10.]
                                neighV[sim][rv][cv]['effS2R'] = updateDictionary(neighV[sim][rv][cv]['effS2R'],tmpNe,tmpDe,numRe,numSe)
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            stTime=1
                            numSs,numRs,tmpDs,tmpNs=[],[],[],[]
                        elif prevS1==4:# prev was transR
                            times[sim][rv][cv]['resIr']+=[intTime/10.]
                            neighV[sim][rv][cv]['resIr'] = updateDictionary(neighV[sim][rv][cv]['resIr'],tmpNi,tmpDi,numRi,numSi)
                            if current==1:#r2S
                                times[sim][rv][cv]['effR2S']+=[eTime/10.]
                                neighV[sim][rv][cv]['effR2S'] = updateDictionary(neighV[sim][rv][cv]['effR2S'],tmpNe,tmpDe,numRe,numSe)
                                eTime=0
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                            stTime=1
                            numSs,numRs,tmpDs,tmpNs=[],[],[],[]

                    else:## in the same state
                         if (current==1) or (current==2):
                             stTime+=1
                             tmpNs+=[updateValue(dN[sim][k],rv,cv,Nrow,Ncol)]
                             tmpDs+=[updateValue(dD[sim][k],rv,cv,Nrow,Ncol)]
                             numSs+=[updateState(states[sim][k],1.,rv,cv,Nrow,Ncol)]
                             numRs+=[updateState(states[sim][k],2.,rv,cv,Nrow,Ncol)]
                         else:#
                             intTime+=1
                             tmpNi+=[updateValue(dN[sim][k],rv,cv,Nrow,Ncol)]
                             tmpDi+=[updateValue(dD[sim][k],rv,cv,Nrow,Ncol)]
                             numSi+=[updateState(states[sim][k],1.,rv,cv,Nrow,Ncol)]
                             numRi+=[updateState(states[sim][k],2.,rv,cv,Nrow,Ncol)]

                    tmpNe+=[updateValue(dN[sim][k],rv,cv,Nrow,Ncol)]
                    tmpDe+=[updateValue(dD[sim][k],rv,cv,Nrow,Ncol)]
                    numSe+=[updateState(states[sim][k],1.,rv,cv,Nrow,Ncol)]
                    numRe+=[updateState(states[sim][k],2.,rv,cv,Nrow,Ncol)]
                    eTime+=1

                    prevS1 = current
                    k+=1
    return [times[sim],neighV[sim],neighS[sim]]

###
def get_results1(ddN,ddD,thresholds):
    ### get the similarity metric
    distD={}
    for sim in ddD:
        datD = ddD[sim]

        distD[sim]=[]

        for k in range(len(datD)):
            distD[sim]+=[getCorr(datD[k])]
    return distD

def get_results2(ddN,ddD,thresholds):
    ### get the average values of properties
    states3,avgR,avgS,avglR,avglS,states2,states={},{},{},{},{},{},{}
    avgER,avgES={},{}

    for sim in ddD:
        datD = ddD[sim]
        datN = ddN[sim]        

        avgR[sim]=[]
        avgS[sim]=[]
        avglR[sim]=[]
        avglS[sim]=[]
        avgER[sim]=[]
        avgES[sim]=[]

        threshN={'S':thresholds['S']['N'],'R':thresholds['R']['N']}
        threshD={'S':thresholds['S']['D'],'R':thresholds['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
        for k in range(len(datD)):
	    tmp1 = states[sim][k]==1
	    tmp2 = states[sim][k]==2
	    tmp3 = states[sim][k]==3
	    tmp4 = states[sim][k]==4
	    tmp5 = states[sim][k]==0
            avgS[sim] +=[np.mean(tmp1)]
            avgR[sim] +=[np.mean(tmp2)]
            avglS[sim] +=[np.mean(tmp3)]
            avglR[sim] +=[np.mean(tmp4)]
            avgES[sim] +=[np.mean((tmp1)*1.+(tmp3)*1.)]
            avgER[sim] +=[np.mean((tmp2)*1.+(tmp4)*1.+(tmp5)*1.)]

    return [avgR,avgS,avglR,avglS,avgER,avgES,states]

def get_results3(ddN,ddD,thresholds):
    ## get the percent of correct contacts
    states2,states={},{}
    states3={}
    contacts={}

    for sim in ddD:
        datD = ddD[sim]
        datN = ddN[sim]        

	contacts[sim]={}
	tmps,tmpo,tmpr=[],[],[]

        threshN={'S':thresholds['S']['N'],'R':thresholds['R']['N']}
        threshD={'S':thresholds['S']['D'],'R':thresholds['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
        for k in range(len(datD)):
            tmp = getCont(states2[sim][k])
            tmps+=[tmp[0]]
            tmpo+=[tmp[1]]
            tmpr+=[tmp[2]]

        contacts[sim]['SS']=tmps
        contacts[sim]['Opp']=tmpo
        contacts[sim]['RR']=tmpr

    return [contacts]

def get_results7(ddN,ddD,thresholds,tfin=-1):
    ## get the escape times
    distD,avgR,avgS,avglR,avglS,states2,states={},{},{},{},{},{},{}
    qres,states3={},{}
    times,neighV,neighS={},{},{}
    for sim in ddD:
        datD = ddD[sim]
        datN = ddN[sim]        

	times[sim]={}
	neighV[sim]={}
	neighS[sim]={}

        threshN={'S':thresholds['S']['N'],'R':thresholds['R']['N']}
        threshD={'S':thresholds['S']['D'],'R':thresholds['R']['D']}
        [states[sim],states2[sim],states3[sim]] = getSimStates(datD,datN,threshD,threshN)    
        [times[sim],neighV[sim],neighS[sim]]=getET_states(states,ddN,ddD,tfin)

    return [times,neighV,neighS]

def getCont(states):

        ss,opp,rr=[],[],[]

	states = np.array(states)
	tmpa= np.roll(states,1,axis=1)#col
	tmpc= np.roll(states,1,axis=0)#row

	## By only doing the two we don't have to correct for over counting
	opp = np.sum(tmpa!=states)+np.sum(tmpc!=states)##+np.sum(tmpd!=states)+np.sum(tmpb!=states)
	ss = np.sum((tmpa==states)*(states==1))+np.sum((tmpc==states)*(states==1))##+np.sum((tmpd!=states)*(states==1))+np.sum((tmpb!=states)*(states==1))
	rr = np.sum((tmpa==states)*(states==2))+np.sum((tmpc==states)*(states==2))##+np.sum((tmpd!=states)*(states==2))+np.sum((tmpb!=states)*(states==2))
	return ss,opp,rr

##############################
def pseudo_potential(N, D, bins):
    #written by Federico Bocci
    '''
    construct the pseudopotential U=-log10(P(N,D))
    from N, D vectors
    '''
    # pipeline to prepare shot noise data for plotting #
    # 1) binning on log scale in [0.1, 10^5]
    # 2) compute 2D histogram of (N, D)
    # 3) change count from 0 to 1 in bins without any count
    # (this is necessary to take the log10 for the pseudopotential)
    # 4) define a 2D matrix 'norm' to account for the uneven binnign when normalizing the PDF
    # 5) to avoid uneven counting of the 'artificial' points with count=1,
    # they are not divided by norm[i][j] but by the maximum possible normalization value norm[-1][-1]

    norm = np.zeros((bins.size - 1, bins.size - 1))
    for i in np.arange(0, bins.size - 1, 1):
        for j in np.arange(0, bins.size - 1, 1):
            norm[i][j] = (bins[i + 1] - bins[i]) * (bins[j + 1] - bins[j])

    dataD = np.ndarray.flatten(D)
    dataN = np.ndarray.flatten(N)
    tmpN=[]
    tmpD=[]
    for i in range(len(dataN)):
        tmpD+=[np.ndarray.flatten(dataD[i])]
        tmpN+=[np.ndarray.flatten(dataN[i])]
    dataN = np.ndarray.flatten(np.array(tmpN))
    dataD = np.ndarray.flatten(np.array(tmpD))
    dist, xbin, ybin = np.histogram2d(dataN, dataD, bins=[bins, bins])

    # add ones to allow computation of U
    dist_new = np.zeros((xbin.size-1, ybin.size-1))
    for i in range(xbin.size - 1):
        for j in range(ybin.size - 1):
            if dist[i][j] == 0:
                dist_new[i][j] = 1.
            else:
                dist_new[i][j] = dist[i][j]

    # normalize probability to 1
    dist_new = dist_new / np.float(np.sum(dist_new))

    # normalize to uneven grid to to get PDF
    for i in range(xbin.size - 1):
        for j in range(ybin.size - 1):
            if dist[i][j] == 0:
                dist_new[i][j] = dist_new[i][j] / norm[-1][-1]
            else:
                dist_new[i][j] = dist_new[i][j] / norm[i][j]
    U = -np.log10(dist_new)
    return U

##############################################################
#
### functions to find minima and __thresholds__ in pseudopotential
#
###

def find_minima(dist, xbins, notch, delta):
    #written by Federico Bocci
    '''
    find the S and R  minima on pseudopotential landscape
    NB this is custom made already knowing that there are only two main minima
    the condition min(notch[i], delta[j])>1 gets rid of adsorbing boundaries around N=0 or D=0
    that can have a high count in some simulations
    '''
    dist_up = np.ones((xbins.size - 1, xbins.size - 1)) * np.amax(dist)
    dist_low = np.ones((xbins.size - 1, xbins.size - 1)) * np.amax(dist)
    for i in range(xbins.size - 1):
        for j in range(xbins.size - 1):
            if i > j and min(notch[i], delta[j])>1.:
                dist_up[i][j] = dist[i][j]
            elif i < j and min(notch[i], delta[j])>1.:
                dist_low[i][j] = dist[i][j]
    i_s, j_s = np.unravel_index(np.argmin(dist_up), dist_up.shape)
    i_r, j_r = np.unravel_index(np.argmin(dist_low), dist_low.shape)
    return i_s, j_s, i_r, j_r


def find_thresholds(x, U, i_s, j_s, i_r, j_r):
    #written by Federico Bocci
    '''
    find distances from S and R minima where the probability is decreased by a 10-fold
    thus the pseudopotential increases by a unit
    '''
    ### find __thresholds__:
    # from Receiver: notch moves left, delta moves up
    Uref = U[i_s][j_s]
    i = i_s - 1
    while (U[i][j_s] - Uref) < 1.:
        i = i - 1
    notch_thr_R = x[i]
    j = j_s + 1
    while (U[i_s][j] - Uref) < 1.:
        j = j + 1
    delta_thr_R = x[j]
    # from Sender: notch moves right, delta moves down
    Uref = U[i_r][j_r]
    i = i_r + 1
    while (U[i][j_r] - Uref) < 1.:
        i = i + 1
    notch_thr_S = x[i]
    j = j_r - 1
    while (U[i_r][j] - Uref) < 1.:
        j = j - 1
    delta_thr_S = x[j]
    return notch_thr_R, delta_thr_R, notch_thr_S, delta_thr_S
