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
    ## based off code written by Federico for finding the minima and getting the threshold at 0.1 of that.

    N= dN[0]
    D= dD[0]

    bins = np.logspace(np.log10(0.1),np.log10(100000),num=121,base=10)
    x = np.zeros(bins.size-1)
    for i in range(x.size):
        x[i] = (bins[i]+bins[i+1])/2.

    U = pseudo_potential(N,D,bins)
    i_s,j_s,i_r,j_r = find_minima(U,bins,x,x)
    ntr,dtr,nts,dts = find_thresholds(x,U,i_s,j_s,i_r,j_r)
    nr = bins[i_s]
    dr = bins[j_s]
    ns = bins[i_r]
    ds = bins[j_r]

    #print('nr,dr,ns,ds')
    #print(nr,dr,ns,ds)
    #print('ntr,dtr,nts,dts')
    #print(ntr,dtr,nts,dts)

    return [nr,dr,ns,ds,ntr,dtr,nts,dts]


#### functions for reading and analyzing data
def convert2timeSeries(dat,size,tstart):
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
    
    [dataN,dataD,dataI] = readFile(filen,tstart)

    myData=[dataN,dataD,dataI]

    return myData

def getCorr(X):
        [nr,nc] = X.shape

        totalC,count=0,0
        for i in range(nr):
                corr,_ = pearsonr(X[i],X[(i+1)%nr])

                corr = (1.-corr)/2 # convert from (-1,1) to (0,1)
                totalC+=corr
                count+=1.
        totalC=totalC/count

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
    ## should be 5 states
    ## 1=Sender, 2=Reciever,3= leftS_state, 4=leftR_state, or 0=none of above    ## go through all the cells in the lattice and find how they behave

    states=[]
    states2=[]
    states3=[]
    iD = tD['R']+ (tD['S']-tD['R'])/2. 
    iN = tN['S']+ (tN['R']-tN['S'])/2. 
    for k in range(len(dD)):
        if k==0:
            ### set the initial as either sender or receiver by breaking the thresholds in halfo
	    tmp = (dD[k]>=dN[k])*1.+(dD[k]<dN[k])*2.
            #tmp = (dD[k]>iD)*(dN[k]<iN)*1. +(dD[k]<iD)*(dN[k]>iN)*2.
        else:
            tmpA = (dD[k]>tD['S'])*(dN[k]<tN['S'])*1. +(dD[k]<tD['R'])*(dN[k]>tN['R'])*2.
            tmp = tmpA + (tmpA==0)*(states[k-1]==1)*3. +(tmpA==0)*(states[k-1]==2)*4. + (tmpA==0)*(states[k-1]==3)*3. + (tmpA==0)*(states[k-1]==4)*4.

        states+=[tmp]
        tmpB = (dD[k]>iD)*(dN[k]<iN)*1. +(dD[k]<iD)*(dN[k]>iN)*2.
        states2+=[tmpB]
        tmpC = (dD[k]>iD)*(dN[k]<iN)*1. -(dD[k]<iD)*(dN[k]>iN)*1.
        states3+=[tmpC]
    return [states,states2,states3]

def getET_states(states,dN,dD,tfin=-1):

    ## first determine how long it is in S or R state
    ## save that in resS or resR
    ## then determine how long in transition state 
    ## if transition completely ot R or S, then 
    ## save Rtime + trans time in etS2R
    ## save Rtime in eS2R
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
                times[sim][rv][cv]={'eS2R':[],'eR2S':[],'bS2R':[],'bR2S':[],'etS2R':[],'etR2S':[],'effS2R':[],'effR2S':[],'resS':[],'resR':[],'resIr':[],'resIs':[]}
                neighV[sim][rv][cv]={'eS2R':{'N':[],'D':[]},'eR2S':{'N':[],'D':[]},'bR2S':{'N':[],'D':[]},'effS2R':{'N':[],'D':[]},'effR2S':{'N':[],'D':[]},'bS2R':{'N':[],'D':[]},'etR2S':{'N':[],'D':[]},'etS2R':{'N':[],'D':[]},'resS':{'N':[],'D':[]},'resR':{'N':[],'D':[]},'resIr':{'N':[],'D':[]},'resIs':{'N':[],'D':[]}}
                neighS[sim][rv][cv]={'eS2R':{'S':[],'R':[]},'eR2S':{'S':[],'R':[]},'effS2R':{'S':[],'R':[]},'effR2S':{'S':[],'R':[]},'bR2S':{'S':[],'R':[]},'bS2R':{'S':[],'R':[]},'etR2S':{'S':[],'R':[]},'etS2R':{'S':[],'R':[]},'resS':{'S':[],'R':[]},'resR':{'S':[],'R':[]},'resIr':{'S':[],'R':[]},'resIs':{'S':[],'R':[]}}
                
		if len(states[sim])>20000:
                	k=10000            
		else:
                	k=1            
                count=0
                prevS1 = states[sim][k-1][rv][cv] ## should be  1-S, 2-R, 3-tS, or 4-tR
                stTime,intTime,eTime=1,1,1
                tmpNs,tmpNi,tmpDs,tmpDi=[],[],[],[]
                numSs,numSi,numRs,numRi=[],[],[],[]
                numSe,numRe,tmpDe,tmpNe=[],[],[],[]

		if tfin!=-1 and len(states[sim])>=tfin:
			final_step =tfin
		else:
			final_step =len(states[sim])
			
                while (k <final_step):
                    current = states[sim][k][rv][cv]
                    if current!=prevS1:
                        if prevS1==1:# prev was S
                            times[sim][rv][cv]['resS']+=[stTime/10.]
                            neighV[sim][rv][cv]['resS']['N']+=[np.mean(tmpNs)]
                            neighV[sim][rv][cv]['resS']['D']+=[np.mean(tmpDs)]
                            neighS[sim][rv][cv]['resS']['R']+=[np.mean(numRs)]
                            neighS[sim][rv][cv]['resS']['S']+=[np.mean(numSs)]
                            if current==2:
				#print 'StoR',stTime,intTime
                                #times[sim][rv][cv]['eS2R']+=[stTime/10.]
                                #neighV[sim][rv][cv]['eS2R']['N']+=[np.mean(tmpNs)]
                                #neighV[sim][rv][cv]['eS2R']['D']+=[np.mean(tmpDs)]
                                #neighS[sim][rv][cv]['eS2R']['R']+=[np.mean(numRs)]
                                #neighS[sim][rv][cv]['eS2R']['S']+=[np.mean(numSs)]
                                times[sim][rv][cv]['effS2R']+=[eTime/10.]
                                neighV[sim][rv][cv]['effS2R']['N']+=[np.mean(tmpNe)]
                                neighV[sim][rv][cv]['effS2R']['D']+=[np.mean(tmpDe)]
                                neighS[sim][rv][cv]['effS2R']['R']+=[np.mean(numRe)]
                                neighS[sim][rv][cv]['effS2R']['S']+=[np.mean(numSe)]
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            intTime=1
                            numSi,numRi,tmpDi,tmpNi=[],[],[],[]
                        elif prevS1==2:# prev was R
                            times[sim][rv][cv]['resR']+=[stTime/10.]
                            neighV[sim][rv][cv]['resR']['N']+=[np.mean(tmpNs)]
                            neighV[sim][rv][cv]['resR']['D']+=[np.mean(tmpDs)]
                            neighS[sim][rv][cv]['resR']['R']+=[np.mean(numRs)]
                            neighS[sim][rv][cv]['resR']['S']+=[np.mean(numSs)]
                            if current==1:
				#print 'RtoS',stTime,intTime
                                #times[sim][rv][cv]['eR2S']+=[stTime/10.]
                                #neighV[sim][rv][cv]['eR2S']['N']+=[np.mean(tmpNs)]
                                #neighV[sim][rv][cv]['eR2S']['D']+=[np.mean(tmpDs)]
                                #neighS[sim][rv][cv]['eR2S']['R']+=[np.mean(numRs)]
                                #neighS[sim][rv][cv]['eR2S']['S']+=[np.mean(numSs)]
                                times[sim][rv][cv]['effR2S']+=[eTime/10.]
                                neighV[sim][rv][cv]['effR2S']['N']+=[np.mean(tmpNe)]
                                neighV[sim][rv][cv]['effR2S']['D']+=[np.mean(tmpDe)]
                                neighS[sim][rv][cv]['effR2S']['R']+=[np.mean(numRe)]
                                neighS[sim][rv][cv]['effR2S']['S']+=[np.mean(numSe)]
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            intTime=1
                            numSi,numRi,tmpDi,tmpNi=[],[],[],[]
                        elif prevS1==3:# prev was transS
                            times[sim][rv][cv]['resIs']+=[intTime/10.]
                            neighV[sim][rv][cv]['resIs']['N']+=[np.mean(tmpNi)]
                            neighV[sim][rv][cv]['resIs']['D']+=[np.mean(tmpDi)]
                            neighS[sim][rv][cv]['resIs']['R']+=[np.mean(numRi)]
                            neighS[sim][rv][cv]['resIs']['S']+=[np.mean(numSi)]
                            if current==2:# s2r
				#print 'IstoR',stTime,intTime
                                #times[sim][rv][cv]['etS2R']+=[intTime/10.]
                                #neighV[sim][rv][cv]['etS2R']['N']+=[np.mean(tmpNi)]
                                #neighV[sim][rv][cv]['etS2R']['D']+=[np.mean(tmpDi)]
                                #neighS[sim][rv][cv]['etS2R']['R']+=[np.mean(numRi)]
                                #neighS[sim][rv][cv]['etS2R']['S']+=[np.mean(numSi)]
                                #times[sim][rv][cv]['eS2R']+= [stTime/10.]
                                #neighV[sim][rv][cv]['eS2R']['N']+=[np.mean(tmpNs)]
                                #neighV[sim][rv][cv]['eS2R']['D']+=[np.mean(tmpDs)]
                                #neighS[sim][rv][cv]['eS2R']['R']+=[np.mean(numRs)]
                                #neighS[sim][rv][cv]['eS2R']['S']+=[np.mean(numSs)]
                                #times[sim][rv][cv]['bS2R']+= [(stTime+intTime)/10.]
                                #neighV[sim][rv][cv]['bS2R']['N']+=[np.mean(np.concatenate([tmpNs,tmpNi]))]
                                #neighV[sim][rv][cv]['bS2R']['D']+=[np.mean(np.concatenate([tmpDs,tmpDi]))]
                                #neighS[sim][rv][cv]['bS2R']['R']+=[np.mean(np.concatenate([numRs,numRi]))]
                                #neighS[sim][rv][cv]['bS2R']['S']+=[np.mean(np.concatenate([numSs,numSi]))]
                                times[sim][rv][cv]['effS2R']+=[eTime/10.]
                                neighV[sim][rv][cv]['effS2R']['N']+=[np.mean(tmpNe)]
                                neighV[sim][rv][cv]['effS2R']['D']+=[np.mean(tmpDe)]
                                neighS[sim][rv][cv]['effS2R']['R']+=[np.mean(numRe)]
                                neighS[sim][rv][cv]['effS2R']['S']+=[np.mean(numSe)]
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                                eTime=0
                            stTime=1
                            numSs,numRs,tmpDs,tmpNs=[],[],[],[]
                        elif prevS1==4:# prev was transR
                            times[sim][rv][cv]['resIr']+=[intTime/10.]
                            neighV[sim][rv][cv]['resIr']['N']+=[np.mean(tmpNi)]
                            neighV[sim][rv][cv]['resIr']['D']+=[np.mean(tmpDi)]
                            neighS[sim][rv][cv]['resIr']['R']+=[np.mean(numRi)]
                            neighS[sim][rv][cv]['resIr']['S']+=[np.mean(numSi)]
                            if current==1:#r2S
				#print 'IrtoS',stTime,intTime
                                #times[sim][rv][cv]['etR2S']+=[intTime/10.]
                                #neighV[sim][rv][cv]['etR2S']['N']+=[np.mean(tmpNi)]
                                #neighV[sim][rv][cv]['etR2S']['D']+=[np.mean(tmpDi)]
                                #neighS[sim][rv][cv]['etR2S']['R']+=[np.mean(numRi)]
                                #neighS[sim][rv][cv]['etR2S']['S']+=[np.mean(numSi)]
                                #times[sim][rv][cv]['eR2S']+= [stTime/10.]
                                #neighV[sim][rv][cv]['eR2S']['N']+=[np.mean(tmpNs)]
                                #neighV[sim][rv][cv]['eR2S']['D']+=[np.mean(tmpDs)]
                                #neighS[sim][rv][cv]['eR2S']['R']+=[np.mean(numRs)]
                                #neighS[sim][rv][cv]['eR2S']['S']+=[np.mean(numSs)]
                                #times[sim][rv][cv]['bR2S']+= [(stTime+intTime)/10.]
                                #neighV[sim][rv][cv]['bR2S']['N']+=[np.mean(np.concatenate([tmpNs,tmpNi]))]
                                #neighV[sim][rv][cv]['bR2S']['D']+=[np.mean(np.concatenate([tmpDs,tmpDi]))]
                                #neighS[sim][rv][cv]['bR2S']['R']+=[np.mean(np.concatenate([numRs,numRi]))]
                                #neighS[sim][rv][cv]['bR2S']['S']+=[np.mean(np.concatenate([numSs,numSi]))]
                                times[sim][rv][cv]['effR2S']+=[eTime/10.]
                                neighV[sim][rv][cv]['effR2S']['N']+=[np.mean(tmpNe)]
                                neighV[sim][rv][cv]['effR2S']['D']+=[np.mean(tmpDe)]
                                neighS[sim][rv][cv]['effR2S']['R']+=[np.mean(numRe)]
                                neighS[sim][rv][cv]['effR2S']['S']+=[np.mean(numSe)]
                                eTime=0
                                numSe,numRe,tmpDe,tmpNe=[],[],[],[]
                            stTime=1
                            numSs,numRs,tmpDs,tmpNs=[],[],[],[]

                    else:## in the same state
                         if (current==1) or (current==2):
                             stTime+=1
                             tmpNs+=[dN[sim][k][(rv-1)%Nrow][cv],dN[sim][k][(rv+1)%Nrow][cv],dN[sim][k][rv][(cv-1)%Ncol],dN[sim][k][rv][(cv+1)%Ncol]]
                             tmpDs+=[dD[sim][k][(rv-1)%Nrow][cv],dD[sim][k][(rv+1)%Nrow][cv],dD[sim][k][rv][(cv-1)%Ncol],dD[sim][k][rv][(cv+1)%Ncol]]
                             numSs+=[(states[sim][k][(rv-1)%Nrow][cv]==1.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==1.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==1.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==1.)*1.]
                             numRs+=[(states[sim][k][(rv-1)%Nrow][cv]==2.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==2.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==2.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==2.)*1.]
                         elif (current==3) or (current==4):
                             intTime+=1
                             tmpNi+=[dN[sim][k][(rv-1)%Nrow][cv],dN[sim][k][(rv+1)%Nrow][cv],dN[sim][k][rv][(cv-1)%Ncol],dN[sim][k][rv][(cv+1)%Ncol]]
                             tmpDi+=[dD[sim][k][(rv-1)%Nrow][cv],dD[sim][k][(rv+1)%Nrow][cv],dD[sim][k][rv][(cv-1)%Ncol],dD[sim][k][rv][(cv+1)%Ncol]]
                             numSi+=[(states[sim][k][(rv-1)%Nrow][cv]==1.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==1.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==1.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==1.)*1.]
                             numRi+=[(states[sim][k][(rv-1)%Nrow][cv]==2.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==2.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==2.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==2.)*1.]

                    tmpNe+=[(dN[sim][k][(rv-1)%Nrow][cv]+dN[sim][k][(rv+1)%Nrow][cv]+dN[sim][k][rv][(cv-1)%Ncol]+dN[sim][k][rv][(cv+1)%Ncol])/4.]
                    tmpDe+=[(dD[sim][k][(rv-1)%Nrow][cv]+dD[sim][k][(rv+1)%Nrow][cv]+dD[sim][k][rv][(cv-1)%Ncol]+dD[sim][k][rv][(cv+1)%Ncol])/4.]
                    numSe+=[(states[sim][k][(rv-1)%Nrow][cv]==1.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==1.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==1.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==1.)*1.]
                    numRe+=[(states[sim][k][(rv-1)%Nrow][cv]==2.)*1.+(states[sim][k][(rv+1)%Nrow][cv]==2.)*1.+(states[sim][k][rv][(cv-1)%Ncol]==2.)*1.+(states[sim][k][rv][(cv+1)%Ncol]==2.)*1.]
                    eTime+=1

		    #print current, prevS1,stTime,intTime
                    prevS1 = current
                    k+=1
    return [times[sim],neighV[sim],neighS[sim]]

def get_results1(ddN,ddD,thresholds):
    distD={}
    for sim in ddD:
        datD = ddD[sim]

        distD[sim]=[]

        for k in range(len(datD)):
            distD[sim]+=[getCorr(datD[k])]
    return distD

def get_results2(ddN,ddD,thresholds):
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
	#tmpb= np.roll(states,-1,axis=1)#col
	tmpc= np.roll(states,1,axis=0)#row
	#tmpd= np.roll(states,-1,axis=0)# row


	## By only doing the two we don't have to correct for over counting
	opp = np.sum(tmpa!=states)+np.sum(tmpc!=states)##+np.sum(tmpd!=states)+np.sum(tmpb!=states)
	ss = np.sum((tmpa==states)*(states==1))+np.sum((tmpc==states)*(states==1))##+np.sum((tmpd!=states)*(states==1))+np.sum((tmpb!=states)*(states==1))
	rr = np.sum((tmpa==states)*(states==2))+np.sum((tmpc==states)*(states==2))##+np.sum((tmpd!=states)*(states==2))+np.sum((tmpb!=states)*(states==2))
	#print opp,ss,rr
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
