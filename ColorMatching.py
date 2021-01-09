# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 18:53:13 2019

@author: idahoo
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 00:22:18 2019

@author: user
"""

import numpy as np
import math
'''
# =============================================================================
#  
# =============================================================================
'''
def XYZ ( E , R , OBS ):
    E=(np.float64(E))
    R=np.float64(R)
    OBS=np.float64(OBS)
    xbar = np.float64(OBS[:,0])
    ybar = np.float64(OBS[:,1])
    zbar = np.float64(OBS[:,2])
    k    = 100 / np.dot(ybar,E)
    X    = k*sum(E*xbar*R)
    Y    = k*sum(E*ybar*R)
    Z    = k*sum(E*zbar*R)
    return(X,Y,Z)
    
def XYZ_New ( E , R , OBS ):
    E=(np.float64(E))
    R=np.float64(R)
    OBS=np.float64(OBS)
    xbar = np.float64(OBS[:,0])
    ybar = np.float64(OBS[:,1])
    zbar = np.float64(OBS[:,2])
    k    = 100 / np.dot(ybar,E)
    X    = (k*sum(E*xbar*R)**(1/3))
    Y    = (k*sum(E*ybar*R)**(1/3))
    Z    = (k*sum(E*zbar*R)**(1/3))
    return(X,Y,Z)
    
    
def DeltaE_New (XYZ1,XYZ2):
    X1      = np.float64(XYZ1[0])
    Y1      = np.float64(XYZ1[1])
    Z1      = np.float64(XYZ1[2])
    X2      = np.float64(XYZ2[0])
    Y2      = np.float64(XYZ2[1])
    Z2      = np.float64(XYZ2[2])
    delta_E  = (((X1-X2)**2)+((Y1-Y2)**2)+((Z1-Z2)**2))**0.5
    return(delta_E)

def LAB_New ( E , R , OBS ):
    E=(np.float64(E))
    R=np.float64(R)
    OBS=np.float64(OBS)
    xbar = np.float64(OBS[:,0])
    ybar = np.float64(OBS[:,1])
    zbar = np.float64(OBS[:,2])
    k    = 100 / np.dot(ybar,E)
    X    = k*sum(E*xbar*(R**(1/3)))
    Y    = k*sum(E*ybar*(R**(1/3)))
    Z    = k*sum(E*zbar*(R**(1/3)))
    XS   = k*sum(E*xbar)
    YS   = 100
    ZS   = k*sum(E*zbar)
    LS   = (116*((Y/YS))**(1/3))-(16)
    aS   = 500*(((X/XS)**(1/3))-((Y/YS))**(1/3))
    bS   = 200*(((Y/YS)**(1/3))-((Z/ZS))**(1/3)) 
    return (LS,aS,bS)

def xy_New(E,R,OBS):
    XYZ1=XYZ_New(E,R,OBS)
    x=XYZ1[0]/sum(XYZ1)
    y=XYZ1[1]/sum(XYZ1)
    return(x,y)
    
def xy(E,R,OBS):
    XYZ1=XYZ(E,R,OBS)
    x=XYZ1[0]/sum(XYZ1)
    y=XYZ1[1]/sum(XYZ1)
    return(x,y)
    
def xy1(X,Y,Z):
    X = np.float64(X)
    Y = np.float64(Y)
    Z = np.float64(Z)
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    return(x,y)    
def LAB ( E , R , OBS ):
    E=(np.float64(E))
    R=np.float64(R)
    OBS=np.float64(OBS)
    xbar = np.float64(OBS[:,0])
    ybar = np.float64(OBS[:,1])
    zbar = np.float64(OBS[:,2])
    k    = 100 / np.dot(ybar,E)
    X    = k*sum(E*xbar*R)
    Y    = k*sum(E*ybar*R)
    Z    = k*sum(E*zbar*R)
    XS   = k*sum(E*xbar)
    YS   = 100
    ZS   = k*sum(E*zbar)
    LS   = (116*((Y/YS))**(1/3))-(16)
    aS   = 500*(((X/XS)**(1/3))-((Y/YS))**(1/3))
    bS   = 200*(((Y/YS)**(1/3))-((Z/ZS))**(1/3)) 
    return (LS,aS,bS)



def LAB1 (E,OBS ,X,Y,Z):
    E=(np.float64(E))
    OBS=np.float64(OBS)
    xbar = np.float64(OBS[:,0])
    ybar = np.float64(OBS[:,1])
    zbar = np.float64(OBS[:,2])
    k    = 100 / np.dot(ybar,E)
    X    = np.float64(X)
    Y    = np.float64(Y)
    Z    = np.float64(Z)
    XS   = k*sum(E*xbar)
    YS   = 100
    ZS   = k*sum(E*zbar)
    LS   = (116*((Y/YS))**(1/3))-(16)
    aS   = 500*(((X/XS)**(1/3))-((Y/YS))**(1/3))
    bS   = 200*(((Y/YS)**(1/3))-((Z/ZS))**(1/3)) 
    return (LS,aS,bS)

def DeltaE (LAB1,LAB2):
    LS1      = np.float64(LAB1[0])
    aS1      = np.float64(LAB1[1])
    bS1      = np.float64(LAB1[2])
    LS2      = np.float64(LAB2[0])
    aS2      = np.float64(LAB2[1])
    bS2      = np.float64(LAB2[2])
    delta_E  = (((LS1-LS2)**2)+((aS1-aS2)**2)+((bS1-bS2)**2))**0.5
    return(delta_E)

'''
# =============================================================================
#  one-constant theory
# =============================================================================
'''
def KOS ( R ):
     R   = np.float64(R)
     KOS = ((1-R)**2)/(2*R)
     return(KOS)
  
def KOSMix(kOs_Sub,kOs_u,C):
    kOs_Sub = np.float64(kOs_Sub)
    kOs_u   = np.float64(kOs_u)
    C       = np.float64(C)
    KOS_Mix = kOs_Sub + np.dot(kOs_u,C)
    return(KOS_Mix)
 
def R(KOS):
    KOS = np.float64(KOS)
    R   = (1+KOS)-(((KOS**2)+2*KOS)**0.5)
    return(R)
'''
# =============================================================================
#  one-constant theory
# =============================================================================
'''
def Wal(L,KOS_M,p,n):
    #enter the length of wavelength as L
    #enter the number of colored sampels as n
    KOS_M   = np.float64(KOS_M)
    p       = np.float64(p)
    KSCOEFS = np.ones([n+1,8])
    KANDS   = np.zeros([L,8])
    OBS     = np.hstack((np.zeros(n),1))
    
    KSCOEFS[:-1,0] = -p.T[:,0]
    KSCOEFS[:-1,1] = -p.T[:,1]
    KSCOEFS[:-1,2] = -p.T[:,2]
    KSCOEFS[:-1,3] = -p.T[:,3]
    
    for i in range(0,16):
        KSCOEFS[:-1,4] = (p[0,:]*KOS_M)[i,:]
        KSCOEFS[:-1,5] = (p[1,:]*KOS_M)[i,:]
        KSCOEFS[:-1,6] = (p[2,:]*KOS_M)[i,:]
        KSCOEFS[:-1,7] = (p[3,:]*KOS_M)[i,:]
        KANDS[i,:]     = np.linalg.pinv(KSCOEFS.T.dot(KSCOEFS)).dot(KSCOEFS.T).dot(OBS)
    return(KANDS)
'''
# =============================================================================
# Spectrophotometric color formulation based on two-constant Kubelka-Munk theory  
# =============================================================================
'''
def S2_C(K,S,STD,L):
    #Calculate the concentration 
    #The *STD symbol is used to summarize instead of KOS_STD 
    #enter the length of wavelength as L
    K     = np.float64(K)
    S     = np.float64(S)
    STD   = np.float64(STD)
    CCOEFS = np.ones([L+1,4])
    OBS   = np.hstack((np.zeros(L),1))
    CCOEFS[:-1,0] = K[:,0] - S[:,0]*STD
    CCOEFS[:-1,1] = K[:,1] - S[:,1]*STD
    CCOEFS[:-1,2] = K[:,2] - S[:,2]*STD
    CCOEFS[:-1,3] = K[:,3] - S[:,3]*STD
    C = np.linalg.pinv(CCOEFS.T.dot(CCOEFS)).dot(CCOEFS.T).dot(OBS)
    return(C)
    
def S2_CW(K,S,STD,L,W):
    #Calculate the concentration Weight matrix
    #The *STD symbol is used to summarize instead of KOS_STD 
    #enter the length of wavelength as L
    W     = np.float64(W)
    K     = np.float64(K)
    S     = np.float64(S)
    STD   = np.float64(STD)
    CCOEFS = np.ones([L+1,4])
    OBS   = np.hstack((np.zeros(L),1))
    CCOEFS[:-1,0] = K[:,0] - S[:,0]*STD
    CCOEFS[:-1,1] = K[:,1] - S[:,1]*STD
    CCOEFS[:-1,2] = K[:,2] - S[:,2]*STD
    CCOEFS[:-1,3] = K[:,3] - S[:,3]*STD
    C = np.linalg.pinv(CCOEFS.T.dot(W).dot(CCOEFS)).dot(CCOEFS.T).dot(W).dot(OBS)
    return(C) 
    
def S2_KAS(C,K,S):
    K  = np.float64(K)
    S  = np.float64(S)
    C  = np.float64(C)
    KOS_Mix = np.float64([sum(i) for i in (C*K)])/np.float64([sum(i) for i in (C*S)])
    sum(C*K)/sum(C*S)
    return(KOS_Mix)
'''
# =============================================================================
# colorimetric color formulation based on two-constant
# =============================================================================
'''
def C2_C(K_STD,R_STD,KANDS,E,OBS):
    #KANDS is for colored sampel
    K_STD  = np.float64(K_STD)
    R_STD  = np.float64(R_STD)
    KANDS  = np.float64(KANDS)
    E      = np.float64(E)
    OBS    = np.float64(OBS)
    S_STD = 1
    K4 = np.vstack((KANDS[:,3]))
    S4 = np.vstack((KANDS[:,7]))
    Phi_K = KANDS[:,0:3]
    Phi_S = KANDS[:,4:7]
    u = np.vstack((1,1,1))
    D_K = np.diag(-2*(R_STD**2)/(1-R_STD**2))
    D_S = np.diag(R_STD*(1-R_STD)/(1+R_STD))
    E=np.diag(E)
    a=np.linalg.pinv(OBS.dot(E).dot(D_K.dot(Phi_K-K4.dot(u.T))+D_S.dot(Phi_S-S4.dot(u.T))))
    C_Cal = a.dot(OBS).dot(E).dot(D_K.dot(K_STD-KANDS[:,3])+D_S.dot(S_STD-KANDS[:,7]))
    return(C_Cal)

import git 
g = git.cmd.Git(git_dir)
g.pull()
git remote add origin https://github.com/IDAHOO75/Color-matching-function.git
git branch -M main
git push -u origin main