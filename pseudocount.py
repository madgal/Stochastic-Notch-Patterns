import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib
from auxFunctions import *

from function_for_si_figures import *

def getContacts_rand(nt,na,tstart,tend):
    contacts=[]
    for sim in range(20):
        fileN = "data/random/Cont_16x16_"+nt+"_n"+str(na)+"_s"+str(sim)+".txt"
        filX=pd.read_csv(fileN).dropna()
        contacts+=list(filX['contOpp'].values[tstart:tend]/5.12)
    return contacts

def getInd(start,stop,by,typeN,binsX,scale=1.):
    xval,yval,zval=[],[],[]
    for i in range(start,stop,by):
    	histX=[]
	results=getContacts_rand(typeN,i,10000,100000)
	for j in range(len(binsX)-1):
		tmp=np.sum((results>=binsX[j])*(results<binsX[j+1]))
		if tmp==0:
			histX+=[1.]
		else:
			histX+=[tmp]
		
	if np.sum(results>=binsX[-1])==0:
		histX+=[1.]
	else:
		histX+=[np.sum(results>=binsX[-1])]

        tmp1,tmp2,tmp3=[],[],[]
        for j in range(len(histX)):
            tmp1+=[i*scale]
            tmp2+=[binsX[j]]
            tmpV = histX[j]*1./np.sum(histX)
            if tmpV!=0:
                tmp3+=[tmpV]#####-np.log10(hist[j]/np.sum(hist))]
            else:
                tmp3+=[1.]
        xval+=[tmp1]
        yval+=[tmp2]
        zval+=[tmp3]
	print i#xval,yval,zval
    
    return xval,yval,zval



'''
binsS = np.arange(50,101,1)
xs,ys,zs= getInd(0,21,1,'shot',binsS)
filen = "dataS12_shot_"
for i in range(len(xs)):
	fileo = open(filen+str(i)+".txt","w")
	fileo.write("x,y,z\n")
	for j in range(len(xs[i])):
		fileo.write("%s,%s,%s\n" %(xs[i][j],ys[i][j],zs[i][j]))
fileo.close()
'''

binsW = np.arange(50,101,1)
xw,yw,zw= getInd(0,210,10,'white',binsW,scale=10)
filen = "dataS12_white_"
for i in range(len(xw)):
	fileo = open(filen+str(i)+".txt","w")
	fileo.write("x,y,z\n")
	for j in range(len(xw[i])):
		fileo.write("%s,%s,%s\n" %(xw[i][j],yw[i][j],zw[i][j]))
fileo.close()
