import re
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def keeperP(name,frame,lines):
    data=open(name,"r")
    dt=data.read().split("\n")
    A=[]
    for line in range(lines*frame+1,lines*(frame+1)):
            D=" ".join(dt[line].split())
            D=D.split(" ")
            D=[float(i) for i in D]
            A.append(D)
            data.close()    
    return np.array(A)

def keeperE(frames,name,lines):
    data=open(name,"r")
    dt=data.read().split("\n")
    B=[]
    for i in range(frames):
        C=" ".join(dt[lines*i].split())
        C=C.split(" ")
        C=[float(i) for i in C]
        B.append(C)
    return np.array(B)

def MaxwellBoltzmann(frame,name,lines):
    B=keeperP(name,frame,lines)
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

def Energy(frames,name,lines):
    sigma=3.4 
    Emm=120.0 
    mass=39.95*1.6747
    D=["K.E","P.E"]
    B=keeperE(frames,name,lines)
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

def rbetween(data):
    A=[]
    B=0
    for i in range(len(data)):
        for j in range(i+1,len(data)):
            B=A[i]-A[j]
            A.append(np.linalg.norm(B))
    return A
        
def RDF(frames,name,lines,stabilization=0):
    A=[]
    for i in range(stabilization,frames):
        A=keeperP(name,stabilization,lines)
    

        
# FRAME WHICH IS IN EQUILIBRIUM
frame=999 # 80/0.005*100
#NUMBER OF LINES FOR EACH FRAME
lines=1001
# NUMBER OF OXYGENS ATOM PER FRAME
name="LA.txt"
#READ DATA
#KEEP THE INFORMATION FROM THIS FRAME
MaxwellBoltzmann(frame,name,lines)
Energy(1000,name,lines)
