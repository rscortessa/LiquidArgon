import re
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def keeperP(dt,frame,lines):
    A=[]
    for line in range(lines*frame+1,lines*(frame+1)):
            D=" ".join(dt[line].split())
            D=D.split(" ")
            D=[float(i) for i in D]
            A.append(D)
            data.close()    
    return np.array(A)

def keeperE(frames,dt,lines):
    B=[]
    for i in range(frames):
        C=" ".join(dt[lines*i].split())
        C=C.split(" ")
        C=[float(i) for i in C]
        B.append(C)
    return np.array(B)

def MaxwellBoltzmann(frame,dt,lines):
    B=keeperP(dt,frame,lines)
    C=B[:,6]
    C=np.sort(C)
    IQR=C[int(3/4*(len(C)-1))]-C[int((len(C)-1)/4)]
    dx=IQR*2/len(C)**(1/3)
    N=int((C[len(C)-1]-C[0])/(dx))+1
    dx=(C[len(C)-1]-C[0])/N
    counts,bins=np.histogram(C,bins=N)
    counts=counts/(len(C)*dx)
    bins=bins-dx/2
    av=np.sum(C)/(len(C))
    var=np.var(C)/np.sqrt(len(C)*(len(C)-1))
    t=av*120/3.0
    print("Mean velocity:",av, "Temperature",t)
    plt.title("Maxwell Boltzmann Distribution"+" $V(T / \epsilon)$")
    plt.ylabel("P(r)")
    plt.xlabel("$r(A^o)$")
    plt.stairs(counts,bins,label="Average velocity="+str(np.round(av,2))+"$\pm$"+str(np.round(var,2)))
    plt.legend()
    plt.savefig("MB.png")

def Energy(frames,dt,lines):
    sigma=3.4 
    Emm=120.0 
    mass=39.95*1.6747
    D=["K.E","P.E"]
    B=keeperE(frames,dt,lines)
    B[:,0]=B[:,0]*sigma*(mass/Emm)**(1/2)
    for KoP in range(2):
        C=B[:,KoP+1]
        plt.figure()
        plt.title(D[KoP]+" $(\epsilon)$")
        plt.ylabel("$Energy(\;\epsilon)$")
        plt.xlabel("$t(10^{-14}s)$")
        if(KoP==0):
            plt.yscale("log")
        plt.plot(B[:,0],C,label=D[KoP])
        plt.legend()
        plt.savefig(D[KoP]+".png")

def PCF(dx,L,dt,lines,frames,stabilization=0):
    mini=0
    maxi=5.0
    N=int((maxi-mini)/dx)
    domine=range(0,frames-stabilization)
    counts=[[] for i in domine]
    for i in domine:
        A=keeperP(dt,i+stabilization,lines)
        A=A[:,0:3]
        B=rbetween(A,L)
        counts[i],bins=np.histogram(B,bins=N,range=(mini,maxi),density=True)
    counts=np.array(counts)
    bins=np.array(bins[1:]-dx/2)
    histogram=[ np.mean(counts[:,i])/(4*np.pi*bins[i]**2*dx) for i in range(N)]
    error=[np.sqrt(np.var(counts[:,i])/((N-1)*(N-2))) for i in range(N)]
    #bins=np.array(bins[1:]-dx/2)
    plt.figure()
    plt.title("Pair Correlation Function "+" $r(\sigma)$")
    plt.ylabel("PCF(r)")
    plt.xlabel("$r(\sigma)$")
    plt.errorbar(bins,histogram,yerr=error,label="Pair Correlation Function")
    plt.axhline(y=1.0,color="r",linestyle="-")
    plt.legend()
    plt.savefig("PCF.png")

    return bins,histogram,error

def rbetween(A,L):
    A=np.array(A)
    m=len(A)
    C=[]
    for i in range(m):
        for j in range(i+1,m):
            D=A[i]-A[j]
            E=np.array([D[i]- round(D[i]/L)*L for i in range(len(D))])
            C.append(np.linalg.norm(E))
    return np.array(C)
# NUMBER OF OXYGENS ATOM PER FRAME
L=10.229
dx=0.002
name="LA.txt"
data=open(name,"r")
dt=data.read().split("\n")         

# FRAME WHICH IS IN EQUILIBRIUM

frame=999 # 80/0.005*100

#NUMBER OF LINES FOR EACH FRAME

lines=1001
#READ DATA
#KEEP THE INFORMATION FROM THIS FRAME

MaxwellBoltzmann(frame,dt,lines)
Energy(1000,dt,lines)
PCF(dx,L,dt,lines,frame,995)
