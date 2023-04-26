import re
import time
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.linear_model import LinearRegression

# Arguments received by console

A=[int(i) for i in sys.argv[1:]]
N,Nwarmup,T,L=A

# Constants
k=1.380649
sigma=3.4 
Emm=120.0 
mass=39.95*1.6747
invdu=sigma*(mass/Emm)**(1/2)
dx=0.01
lines=1001
stepdecay=50
dv=(k*Emm/mass)**(1/2)
deltat=0.01
# KeeperP collects the data of a file.

# Before using the function it is required to have created a list containing the lines of the file
# That is possible following the next commands:
# data=open(name,"r") name is the name of the file
# dt=data.read().split("\n") dt is the list which its elements are the lines of the file
# "start" refers to the inital line where the data is begin to be extracted
# "amount" refers to the amount of lines to be extracted
# step is the sequence of the extraction (The extraction happens every "step" lines) 
# (mini,maxi) define the column range extracted from the file

def keeperP(dt,start,amount,step,mini,maxi):
    A=[]
    for line in range(start,start+amount,step):
            D=" ".join(dt[line].split())
            D=D.split(" ")
            D=[float(i) for i in D[mini:maxi]]
            A.append(D)
            data.close()    
    return np.array(A)

# MaxwellBoltzmann computes the probability distribution of the magnitude of the velocity |V|

# frame is the time step in which the distribution is computed
# lines is the number of lines printed in a time step (N+1) where N is the number of particles 
# dt is the list which contains all the lines of the file
# The Diaconis rule is taking into account to generate the histogram

def MaxwellBoltzmann(frame,dt,lines):

    B=keeperP(dt,frame*lines+1,lines,1,7,8)
    B=np.array(B[:,0])
    C=np.sort(B)
    C=np.array(C)

    #Diaconis rule
    
    IQR=C[int(3/4*(len(C)-1))]-C[int((len(C)-1)/4)]
    dx=IQR*2/len(C)**(1/3)
    N=int((C[len(C)-1]-C[0])/(dx))+1
    dx=(C[len(C)-1]-C[0])/N

    #Create the histogram

    counts,bins=np.histogram(C,bins=N)
    counts=counts/(len(C)*dx)
    bins=bins-dx/2
    av=np.sum(C)/(len(C))

    #Calculate the error

    var=np.var(C)/np.sqrt(len(C)*(len(C)-1))

    #Calculate the velocity multiplying by the right constants

    t=av*120/3.0
    print("Mean velocity:",av, "Temperature",t," K")
    plt.figure()
    plt.title("Maxwell Boltzmann Distribution"+" $V [(\epsilon / M)^{1/2}]$ \n "+"$L="+str(L/10.0)+"\sigma \; T="+str(T/10.0)+"K\; N_{warmup}="+str(Nwarmup)+"$")
    plt.ylabel("$P(|V|)$")
    plt.xlabel("$V [(\epsilon / M)^{1/2}]$")
    plt.stairs(counts,bins,label="Average velocity="+str(np.round(av,2))+"$\pm$"+str(np.round(var,3)))
    #plt.plot(x,y,label="Maxwell Boltzmann "+str(T/10.0))
    plt.legend()
    plt.savefig("MB.png")

# Energy. The function plots the potential and kinetic energy as a function of the time step.

# frames is the number of time steps of the simulation
# lines is the number of lines printed in a time step (N+1) where N is the number of particles 
# dt is the list which contains all the lines of the file
# du is the normalization of the time

    
def Energy(frames,dt,lines,du,i,step):
    D=["Kinetic Energy","Potential Energy"]
    B=keeperP(dt,0,frames*step,step,0,3)
    B[:,0]=B[:,0]*du

    # if KoP is 0, the kinetic energy is plotted, in the other hand with KoP=1 is the potential energy  
    for KoP in range(2):
        
        C=B[:,KoP+1]
        plt.figure()
        plt.title(D[KoP]+"$[\epsilon]$ \n "+"$L="+str(L/10.0)+"\sigma \; T="+str(T/10.0)+"K\; N_{warmup}="+str(Nwarmup)+"$")
        plt.ylabel("$Energy(\;\epsilon)$")
        plt.xlabel("$t(10^{-12}s)$")
        
        if(KoP==0):
            plt.yscale("log")
        plt.plot(B[:,0],C,label=D[KoP])
        plt.legend()
        plt.tight_layout()
        plt.savefig(i+D[KoP]+".png")

#PCF1 evaluates the pair correlation function.
# dx is the width of the discretization of the space (Used to compute the density of particles n(r))
# L is the length of the system simulated
# Nu is the number of the particles implemented in the simulation.
# dt is the list which contains all the lines of the file.
# "lines" is the number of lines of the file for each sample in total they are Nu*(Nu-1)/2 (number of links between particles).
# frames is the number of samples.

def PCF1(dx,L,Nu,dt,lines,frames):

    #Initialize the variables
    #Store the data of the relative positions in A 
    L=L/(10.0)
    mini=0
    maxi=L/2
    N=int((maxi-mini)/dx)
    A=np.array(keeperP(dt,0,frames*lines,1,0,1))
    iteration=range(frames)
    counts=[[] for i in iteration]
    bins=[]

    #Create a histogram for each sample, also the pair correlation function is calculated

    for i in iteration:
        counts[i],bins=np.histogram(A[lines*i:lines*(i+1)],bins=N,range=(mini,maxi),density=False)
        counts[i]=[ 2*counts[i][j]/(4*np.pi*(bins[j]+dx/2.0)**2*dx*(Nu**2/L**3)) for j in range(N)]

    # It is performed an average and the error is obtained    

    counts=np.array(counts)
    RDF=[np.mean(counts[:,i]) for i in range(N)]
    bins=np.array(bins[1:]-dx/2)
    error=[np.sqrt(np.var(counts[:,i])/(frames-1)) for i in range(N)]

    # The data is plotted

    plt.figure()
    plt.title("$Pair\; Correlation\; Function$ "+" $G(r)$ \n "+"$L="+str(L)+"\sigma \; T="+str(T/10.0)+"K\; N_{warmup}="+str(Nwarmup)+"$")
    plt.ylabel("$G(r)$")
    plt.xlabel("$r/ \sigma$")
    plt.errorbar(bins,RDF,yerr=error,xerr=dx,fmt=".k",ecolor="blue",label="$G(r)$",ms=1)
    plt.axhline(y=1.0,color="r",linestyle="-")
    plt.legend()
    plt.savefig("PCF.pdf")

# It is calculated the autocorrelation of the velocity:
# frame is the number 

def AutocorrelationV(frame,dt,lines,rep,timespace,stabilization=0):
    
    D=range(0,frame-stabilization)
    A=range(N)
    t=[i for i in D]
    B=[[] for i in A]
    difussion=[0,0]
    titles="Velocity Autocorrelation $ < V(0)\cdot V(t) > $\n "+"$L="+str(L/10.0)+"\sigma \; T="+str(T/10.0)+"K\; N_{warmup}="+str(Nwarmup)+"$"
      
    ylabels="$Velocity\; Autocorrelation \; [V(0)^2]$"
    save="AutocorrelationV.pdf"
    ra2=[[] for h in range(rep)]

    for h in range(rep):
        C=keeperP(dt,(stabilization+h*timespace)*lines+1,lines*(frame-stabilization),lines,3,6)            
        for i in A:
            #E=[ np.dot(np.array(C[i,0:3]),np.array(C[0,0:3]))/(np.dot(np.array(C[0,0:3]),np.array(C[0,0:3]))) for i in D]
            E=[ np.dot(np.array(C[i,0:3]),np.array(C[0,0:3])) for i in D]
            B[i]=E
            
        B=np.array(B)
        ra2[h]=[np.mean(B[:,k]) for k in D]

    ra2=np.array(ra2)    
    r2=[np.mean(ra2[:,k])/(np.mean(ra2[:,0])) for k in D]
    error=[np.sqrt(np.var(ra2[:,k])/(rep-1))/np.mean(ra2[:,0]) for k in D]
    
    for x in range(len(r2)):
        difussion[0]=difussion[0]+r2[x]
        difussion[1]=difussion[1]+error[x]

    difussion[0]=difussion[0]/3*(dv)**2*(deltat)*(T/Emm) #10^{-8}m^2s^{-1} --> 10^{-4}cm^2s^{-1} 
    difussion[1]=difussion[1]/3*(dv)**2*(deltat)*(T/Emm)
    print("The difussion coefficient is:",difussion[0],"$\pm$",difussion[1])
    
    plt.figure()
    plt.title(titles)
    plt.ylabel(ylabels)
    #plt.text(25,0.8,"$D_0=$"+str(round(difussion[0],3))+"$\pm$"+str(round(difussion[1],3))+" $[10^{-4}cm^2s^{-1}]$")
    plt.xlabel("$ t(10^{-12}s)$")
    plt.errorbar(t,r2,yerr=error,fmt=".k",ecolor="blue",label=ylabels)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save)



name="WARMUPCHECK-N"+str(N)+"Nw"+str(Nwarmup)+"T"+str(T)+"L"+str(L)+".txt"
data=open(name,"r")
dt=data.read().split("\n")

Energy(Nwarmup,dt,lines,invdu,"W",1)

name="HISTORY-N"+str(N)+"Nw"+str(Nwarmup)+"T"+str(T)+"L"+str(L)+".txt"
data=open(name,"r")
dt=data.read().split("\n")

MaxwellBoltzmann(0,dt,lines-1)
AutocorrelationV(150,dt,lines,N,stepdecay,0)
Energy(N*stepdecay,dt,lines,invdu,"M",lines)

name="RPOSITIONS-N"+str(N)+"Nw"+str(Nwarmup)+"T"+str(T)+"L"+str(L)+".txt"
data=open(name,"r")
dt=data.read().split("\n")
Npart=1000
lines=int(Npart/2*(Npart-1))
if(N>10):
    N=10
PCF1(dx,L,Npart,dt,lines,N)
