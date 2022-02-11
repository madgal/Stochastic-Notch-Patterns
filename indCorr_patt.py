import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib as mpl
from scipy.optimize import curve_fit

mpl.rcParams.update({'font.size':20})
def func(x,a,b):
        return a*np.exp(-x*b)
def func2(x,a,b,c):
        return a*x**b+c
def func3(x,a,b,a0,b0,c):
        return a*np.exp(-x*b)+a0*np.exp(-x*b0)
def func4(x,a,b,c):
        return a*np.exp(-x*b)
def func5(x,a,b,c):
        return a*np.exp(-x*b)
def func6(x,a,b,c):
        return a*x**b+c

def leftICS(fil2,tstart,tf1):

    x= np.arange(tstart,tf1)/10.
    y = np.array(fil2['qres'][tstart:tf1])

    try:
        popt3, pcov = curve_fit(func3,x,y)
    except:
        popt3=[]
    try:
        popt4, pcov = curve_fit(func4,x[:2500],y[:2500])
    except:
        popt4=[]
    try:
        popt5, pcov = curve_fit(func6,x[:2500],y[:2500])
    except:
        popt5=[]

    return popt3,popt4,popt5

def arriveFinal(filX,tstart):
    r2 = np.max(filX['qres2'][tstart:])- np.min(filX['qres2'][tstart:])
    finalV = np.abs(filX['qres2']-filX['qres2'][len(filX['qres2'])-1])

    x= np.arange(tstart,len(filX['qres2']))/10.
    y = np.array(filX['qres2'][tstart:])

    try:
        popt3, pcov = curve_fit(func3,x,y)
    except:
        popt3=[]
    try:
        popt4, pcov = curve_fit(func6,x[-50000:],y[-50000:])
    except:
        popt4=[]
    try:
        popt5, pcov = curve_fit(func5,x[1000:],y[1000:])
    except:
        popt5=[]


    return popt3,popt5,popt4
def setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na):
        ax1.set_title("Compared with initial pattern")
        ax2.set_title("Compared with final pattern")

        plt.suptitle("$\sigma_{"+nt+"}$="+str(na)+"\n",fontsize=30)
        ax1.set_xticks([])
        ax2.set_xticks([])
        ax4.legend(loc='lower right')
        ax5.legend(loc='lower right')
        ax7.legend(loc='lower right')
        ax8.legend(loc='lower right')
        ax7.set_xlabel("Time (hr)")
        ax8.set_xlabel("Time (hr)")
        ax1.set_ylabel("q")
        ax2.set_ylabel("q")
        ax4.set_ylabel("q")
        ax5.set_ylabel("q")
        ax7.set_ylabel("q")
        ax8.set_ylabel("q")


########
########
########
def p0_s0():
    nt,na,ext='shot',0,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

################
################
################
def p0_s5():
    nt,na,ext='shot',5,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p0_s10():
    nt,na,ext='shot',10,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p0_s15():
    nt,na,ext='shot',15,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p0_s20():
    nt,na,ext='shot',20,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p0_w0():
    nt,na,ext='white',0,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p0_w50():
    nt,na,ext='white',50,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p0_w100():
    nt,na,ext='white',100,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########

def p0_w150():
    nt,na,ext='white',150,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p0_w200():
    nt,na,ext='white',200,'p0'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

##### D1500########
########
########
def p4_s0():
    nt,na,ext='shot',0,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

################
################
################
def p4_s5():
    nt,na,ext='shot',5,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p4_s10():
    nt,na,ext='shot',10,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p4_s15():
    nt,na,ext='shot',15,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p4_s20():
    nt,na,ext='shot',20,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p4_w0():
    nt,na,ext='white',0,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p4_w50():
    nt,na,ext='white',50,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p4_w100():
    nt,na,ext='white',100,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########

def p4_w150():
    nt,na,ext='white',150,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p4_w200():
    nt,na,ext='white',200,'p4'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#### D5000 ####
def p7_s0():
    nt,na,ext='shot',0,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

################
################
################
def p7_s5():
    nt,na,ext='shot',5,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p7_s10():
    nt,na,ext='shot',10,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p7_s15():
    nt,na,ext='shot',15,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p7_s20():
    nt,na,ext='shot',20,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p7_w0():
    nt,na,ext='white',0,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p7_w50():
    nt,na,ext='white',50,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p7_w100():
    nt,na,ext='white',100,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########

def p7_w150():
    nt,na,ext='white',150,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p7_w200():
    nt,na,ext='white',200,'p7'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

### p8###
def p8_s0():
    nt,na,ext='shot',0,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

################
################
################
def p8_s5():
    nt,na,ext='shot',5,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p8_s10():
    nt,na,ext='shot',10,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p8_s15():
    nt,na,ext='shot',15,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p8_s20():
    nt,na,ext='shot',20,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p8_w0():
    nt,na,ext='white',0,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

#########
#########
#########
def p8_w50():
    nt,na,ext='white',50,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p8_w100():
    nt,na,ext='white',100,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1])-0.001,labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1])-0.001,labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:])-0.001,labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:])-0.001,labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########

def p8_w150():
    nt,na,ext='white',150,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
#########
#########
#########
def p8_w200():
    nt,na,ext='white',200,'p8'
    Dir='patt'
    tstart,tf1=1,10000

    fileN = "figures_mar2021/"+Dir+"/data/ferr2_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    fileN2 = "figures_mar2021/"+Dir+"/data/ferr_16x16_"+str(nt)+"_n"+str(na)+"_"+str(ext)+"_s0.txt"
    roundVal=6
    if os.path.exists(fileN) and os.path.exists(fileN2):
        filX=pd.read_csv(fileN).dropna()
        fil2=pd.read_csv(fileN2).dropna()

        p3,p4,p6 = leftICS(fil2,tstart,tf1)
        pf3,pf5,pf6 = arriveFinal(filX,tstart)

        fig,((ax1,ax2),(ax4,ax5),(ax7,ax8)) = plt.subplots(3,2,figsize=(30,15))
        fig.subplots_adjust(hspace=0.01,wspace=0.05)
        xv1 = np.arange(0,len(fil2['qres'][tstart:tf1]))*0.1
        xv = np.arange(0,len(filX['qres'][tstart:]))*0.1

        ax1.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        ax4.plot(xv1,fil2['qres'][tstart:tf1],linewidth=4,label='data')
        if len(p3)>0:
                labn=str(round(p3[0],2))+"exp(-"+str(round(p3[1],roundVal))+"x)+"+str(round(p3[2],2))+"exp(-"+str(round(p3[3],roundVal))+"x)"
                ax4.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax4.plot(xv1,func3(xv1,*p3),linewidth=4,label='exp')
		tmp1=p3[1]
		tmp2=p3[3]
	else:
		tmp1=np.nan
		tmp2=np.nan
        ax7.plot(xv1,fil2['qres'][tstart:tf1],linewidth=6,label='data')
        if len(p6)>0:
                labn=str(round(p6[0],roundVal))+"x^"+str(round(p6[1],roundVal))+"+"+str(round(p6[2],2))
                ax7.text(10,np.max(fil2['qres'][tstart:tf1]),labn)
                ax7.plot(xv1,func6(xv1,*p6),linewidth=2,label='pst')
		tmp3=p6[1]
	else:
		tmp3=np.nan

        ax2.plot(xv,filX['qres2'][tstart:])
        ax5.plot(xv,filX['qres2'][tstart:],label='data')
        if len(pf3)>0:
                labn=str(round(pf3[0],2))+"exp(-"+str(round(pf3[1],roundVal))+"x)+"+str(round(pf3[2],2))+"exp(-"+str(round(pf3[3],roundVal))+"x)"
                ax5.text(100,np.max(filX['qres2'][tstart:]),labn)
                ax5.plot(xv,func3(xv,*pf3),linewidth=4,label='exp')
		tmp4=pf3[1]
		tmp5=pf3[3]
	else:
		tmp4=np.nan
		tmp5=np.nan
        ax8.plot(xv,filX['qres2'][tstart:],label='data',linewidth=6)
        if len(pf6)>0:
                labn=str(round(pf6[0],roundVal))+"x^"+str(round(pf6[1],roundVal))+"+"+str(round(pf6[2],2))
                ax8.text(10,np.max(filX['qres2'][tstart:]),labn)
                ax8.plot(xv[1000:],func6(xv[1000:],*pf6),linewidth=2,label='pft')
		tmp6=pf6[1]
	else:
		tmp6=np.nan
                
                
        setL(ax1,ax2,ax4,ax5,ax7,ax8,nt,na)

        plt.show()
        fig.savefig("corrTimes_"+str(nt)+"_n"+str(na)+"_"+str(ext)+".png",bbox_inches='tight')
	return tmp1,tmp2,tmp3,tmp4,tmp5,tmp6


