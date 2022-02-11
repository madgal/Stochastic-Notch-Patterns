import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from basicFunctions import *
import time
from mpl_toolkits import mplot3d

def getContacts_rand(nt,na,tstart):
    contacts=[]
    fileN = "~/Research/signaling/paper1/figures_Oct2020/data/Cont_16x16_"+nt+"_n"+str(na)+"_s0.txt"
    filX=pd.read_csv(fileN).dropna()
    return filX['contOpp'].values[tstart:]

xval,yval,zval=[],[],[]
for i in range (20):
    hist,bins = np.histogram(getContacts_rand('shot',i,0),bins=10)
    #print hist,bins
    tmp1,tmp2,tmp3=[],[],[]
    for j in range(len(hist)):
        tmp1+=[i]
        tmp2+=[bins[j]]
        tmp3+=[hist[j]*1./np.sum(hist)]#-np.log10(hist[j]/np.sum(hist))]
    xval+=[tmp1]
    yval+=[tmp2]    
    zval+=[tmp3]    

fig = plt.figure()
ax = plt.axes(projection='3d')
zval = -np.log10(zval)
ax.contour3D(xval, yval, zval, 50, cmap='Reds_r')
ax.set_xlabel('$\sigma_{shot}$')
ax.set_ylabel('Contacts')
ax.set_zlabel('Probability')
ax.set_title('Random Initial Lattice')
plt.show()

fig = plt.figure()
zval = -np.log10(zval)
plt.contourf(xval, yval, zval,50, cmap='rainbow')
plt.xlabel('$\sigma_{shot}$')
plt.ylabel('Contacts')
plt.show()

xval,yval,zval=[],[],[]
for i in range (20):
    hist,bins = np.histogram(getContacts_rand('shot',i,10000),bins=10)
    #print hist,bins
    tmp1,tmp2,tmp3=[],[],[]
    for j in range(len(hist)):
        tmp1+=[i]
        tmp2+=[bins[j]]
        tmp3+=[hist[j]*1./np.sum(hist)]#-np.log10(hist[j]/np.sum(hist))]
    xval+=[tmp1]
    yval+=[tmp2]    
    zval+=[tmp3]    

fig = plt.figure()
ax = plt.axes(projection='3d')
zval = -np.log10(zval)
ax.contour3D(xval, yval, zval,50,  cmap='Reds_r')
ax.set_xlabel('$\sigma_{shot}$')
ax.set_ylabel('Contacts')
ax.set_zlabel('Probability')
ax.set_title('Random Initial Lattice-after relaxation')
plt.show()
plt.close()

fig = plt.figure()
zval = -np.log10(zval)
plt.contourf(xval, yval, zval,50, cmap='rainbow')
plt.colorbar(label='Pseudopotential -log10(P(contacts))')
plt.xlabel('$\sigma_{shot}$')
plt.ylabel('Contacts')
plt.title('Random Initial Lattice-after relaxation')
plt.show()
exit()




def getContacts_check(nt,na,tstart):
    contacts=[]
    fileN = "~/Research/signaling/paper1/figures_mar2021/patt/data/Cont_16x16_"+nt+"_n"+str(na)+"_p8_s0.txt"
    filX=pd.read_csv(fileN).dropna()
    return filX['contOpp'].values[tstart:]


### checkerboard
xval,yval,zval=[],[],[]
for i in range (20):
    hist,bins = np.histogram(getContacts_check('shot',i,0),bins=10)
    #print hist,bins
    tmp1,tmp2,tmp3=[],[],[]
    for j in range(len(hist)):
        tmp1+=[i]
        tmp2+=[bins[j]]
        tmp3+=[hist[j]*1./np.sum(hist)]#-np.log10(hist[j]/np.sum(hist))]
    xval+=[tmp1]
    yval+=[tmp2]    
    zval+=[tmp3]    

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(xval, yval, zval, 50, cmap='binary')
ax.set_xlabel('$\sigma_{shot}$')
ax.set_ylabel('Contacts')
ax.set_zlabel('Probability')
ax.set_title('Checkerboard Initial Lattice')
plt.show()






def getContacts_nuc(nt,na,tstart):
    contacts=[]
    fileN = "~/Research/signaling/paper1/figures_mar2021/patt/data/Cont_16x16_"+nt+"_n"+str(na)+"_p7_s0.txt"
    filX=pd.read_csv(fileN).dropna()
    return filX['contOpp'].values[tstart:]

xval,yval,zval=[],[],[]
for i in range (20):
    try:
        hist,bins = np.histogram(getContacts_nuc('shot',i,0),bins=10)
        tmp1,tmp2,tmp3=[],[],[]
        for j in range(len(hist)):
            tmp1+=[i]
            tmp2+=[bins[j]]
            tmp3+=[hist[j]*1./np.sum(hist)]#-np.log10(hist[j]/np.sum(hist))]
        xval+=[tmp1]
        yval+=[tmp2]    
        zval+=[tmp3]    
    except:
	print "error for shot="+str(i)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(xval, yval, zval, 50, cmap='binary')
ax.set_xlabel('$\sigma_{shot}$')
ax.set_ylabel('Contacts')
ax.set_zlabel('Probability')
ax.set_title('Nucleating Initial Lattice')
plt.show()
