# -*- coding: utf-8 -*-
"""

@author: tnm12

Simulator 2 version 0.1.6

Patch Notes 2.0.1.0:
    1.) Converted rate constants from scalars to N sized vectors / NxN matrices
    2.) Added global rate constants function to change rate constants
    3.) Added capability to change individual rate constants to molecular species
    
Patch Notes 2.0.1.1:
    1.) Added in reverse rate constant for the reporter (krepr)
    2.) Using this new rate constant, added in reporter reverse reaction equations to dR, dS, dO, and dRO
    4.) Changed rate constant if statements from False to 'False' so that a constant could be made 0
    5.) Created krsdGmcsd,krevOmcsd variables so rate equations could be correctly updated to handle varying rate constants across rows
    
Patch Notes 2.0.1.2:
    1.) Added threshold reactions
    
    
Patch Notes 2.0.1.3:
    1.) Correctly added threshold reactions
    2.) Converted krev terms along the diagonal to 0 making the equations more accurate 
    3.) Allows user to specify sovler_ivp method
    
Patch Notes 2.0.1.4:
    1.) Added fuel reactions
    2.) Found rounding issue, so took away dtype=int from initializations
    
Patch Notes 2.0.1.5:
    1.) Added AND gate reactions
    2.) Organized Rate Equations
    
Patch Notes 2.0.1.6:
    1.) Added comparator gate reactions
    2.) Original simulator had two incorrect terms in WTA equation
    
Patch Notes 2.0.1.7:
    1.) Attempted to switch to numbakit-ode
    
Patch Notes 2.0.1.8:
    1.) Added in different transcripton matrices, and ability to change these transcription rates

Patch Notes 2.0.1.9:
    1.) Added degradation reactions
    
Patch Notes 2.0.2.0:
    1.) Added fuel reactions for AND gates

Patch Notes 2.0.2.1:
    1.) Added discontinuous function
    2.) Added transcription calibration function


"""

import numpy as np
import scipy.integrate as spi
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import re
import time
import math


def rate_eqs(t,y,ktxnO,ktxnG,ktxnTh,ktxnF,ktxnAG,ktxnCG,krz,krsd,krev,krep,krepr,kth,krzTh,krsdF,krevF,krevA,krsdA,krzA,krevCG,krsdCG,krzCG,leak,leakA,Otempm,Gtempm,Thtempv,Ftempv,AGtempm,AGmap,CGtempm,CGmap,kssdO,kssdF,kdsduG,kdsdG,kdsdGO,kdsduAG,kdsdAG,kdsduCG,kdsdCG,kdrd,N):

    
    uG = y[0:N**2]
    G = y[N**2:2*N**2]
    O = y[2*N**2:3*N**2]
    GO = y[3*N**2:4*N**2]
    Rv = y[4*N**2:N+4*N**2]
    Sv = y[N+4*N**2:2*N+4*N**2]
    RO = y[2*N+4*N**2:2*N+5*N**2]
    uThv = y[2*N+5*N**2:3*N+5*N**2]
    Thv = y[3*N+5*N**2:4*N+5*N**2]
    Fv = y[4*N+5*N**2:5*N+5*N**2]
    GFv = y[5*N+5*N**2:6*N+5*N**2]
    uAG = y[6*N+5*N**2:6*N+6*N**2]
    AG = y[6*N+6*N**2:6*N+7*N**2]
    AGOa = y[6*N+7*N**2:6*N+8*N**2]
    AGOb = y[6*N+8*N**2:6*N+9*N**2]
    uCG = y[6*N+9*N**2:6*N+10*N**2]
    CG = y[6*N+10*N**2:6*N+11*N**2]
    CGOa = y[6*N+11*N**2:6*N+12*N**2]
    CGOb = y[6*N+12*N**2:6*N+13*N**2]
    AGFbv = y[6*N+13*N**2:7*N+13*N**2]
    

    ##Reshaping
    
    uGm = np.array(uG).reshape(N,N)
    Gm = np.array(G).reshape(N,N)
    Om = np.array(O).reshape(N,N)
    
    
    GOm = np.array(GO).reshape(N,N)
    
    Rm = np.diag(Rv)
    
    Sm = np.diag(Sv)
    
    ROm = np.array(RO).reshape(N,N)
    
    Thm = np.diag(Thv)
    
    GFm = np.diag(GFv)
    
    uAGm = np.array(uAG).reshape(N,N)
    AGm = np.array(AG).reshape(N,N)
    AGOam = np.array(AGOa).reshape(N,N)
    AGObm = np.array(AGOb).reshape(N,N)
    
    AGFbm = np.diag(AGFbv)
    
    uCGm = np.array(uCG).reshape(N,N)
    CGm = np.array(CG).reshape(N,N)
    CGOam = np.array(CGOa).reshape(N,N)
    CGObm = np.array(CGOb).reshape(N,N)
    
    
    
    ##Sum Vectors/Matrices
    
    krsd_Gmcsv = np.sum(krsd*Gm,axis=1)
    krsd_Gmcsd = np.diag(krsd_Gmcsv)
    
    krev_Omcsv = np.sum(krev*Om,axis=1)
    krev_Omcsd = np.diag(krev_Omcsv)
    
    krevF_Omrsv = np.sum(krevF*Om,axis=0)
    
    Omrsv = np.sum(Om,axis=0)
    Omrsd = np.diag(Omrsv)
    
    GOmrsv = np.sum(GOm,axis=0)
    GOmrsd = np.diag(GOmrsv)
    
    ROmrsv = np.sum(ROm,axis=0)
    
    krsdFm = np.diag(krsdF*Fv)
    
    krsdA_AGmcsv = np.sum(AGm*krsdA,axis=1)
    krsdA_AGmcsd = np.diag(krsdA_AGmcsv)
    
    krsdA_AGOamcsv = np.sum(AGOam*krsdA,axis=1)
    krsdA_AGOamcsd = np.diag(krsdA_AGOamcsv)
    
    krevA_Omcsv = np.sum(Om*krevA,axis=1)
    krevA_Omcsd = np.diag(krevA_Omcsv)
    
    AGObmrsv = np.sum(AGObm,axis=0)
    AGObmrsd = np.diag(AGObmrsv)
    
    krsdCG_CGmrsv = np.sum(CGm*krsdCG,axis=0)
    krsdCG_CGmrsd = np.diag(krsdCG_CGmrsv)
    
    krsdCG_CGmcsv = np.sum(CGm*krsdCG,axis=1)
    krsdCG_CGmcsd = np.diag(krsdCG_CGmcsv)
    
    krevCG_CGOamrsv = np.sum(CGOam*krevCG,axis=0)
    krevCG_CGOamrsd = np.diag(krevCG_CGOamrsv)
    
    krevCG_CGObmrsv = np.sum(CGObm*krevCG,axis=0)
    krevCG_CGObmrsd = np.diag(krevCG_CGObmrsv)
    
    CGannA = np.diag(np.sum((Om @ CGmap.T),axis=0))
    CGannB = np.diag(np.sum((Om @ CGmap),axis=0))
    
    OMannA = Om @ np.diag(np.sum(krsdCG*CGOam @ CGmap,axis=0))
    OMannB = Om @ np.diag(np.sum(krsdCG*CGObm @ CGmap.T,axis=0))
    
    kdrd_ROmrsv = np.sum(ROm*kdrd,axis=0)
    
    ##Rate Equations
    duGm = (-krz*uGm + ktxnG*Gtempm - kdsduG*uGm).flatten()
    dGm = (krz*uGm - (Omrsd @ (krsd*Gm))  + GOmrsd @ (krev*Om) - kdsdG*Gm).flatten()
    dOm = (-kssdO*Om\
           -Om @ krsdCG_CGmrsd + krevCG*CGObm - Om @ krsdCG_CGmcsd + krevCG*CGOam - OMannB - OMannA \
           + Omrsd @ (krsdA*AGOam) -Om @ krsdA_AGmcsd - Om @ krsdA_AGOamcsd  + AGObmrsd @ (-krevA*Om) + AGObm @ (krevA_Omcsd) + AGmap @ (leakA*ktxnAG*AGtempm) \
           + GOm @ krsdFm - (krevF*Om) @ GFm + AGObm @ krsdFm - (krevF*Om) @ AGFbm\
           -kth*Om @ Thm \
           + ROm @ (krepr*Sm) + Omrsd @ (krsd*Gm) - Om @ krsd_Gmcsd - Om @ (krep*Rm) + GOmrsd @ (-krev*Om) + GOm @ krev_Omcsd  + leak*ktxnG*Gtempm + ktxnO*Otempm).flatten()
    dGOm = (-kdsdGO*GOm \
            -GOm @ krsdFm + (krevF*Om) @ GFm \
            + Om @ krsd_Gmcsd - (GOm) @ krev_Omcsd ).flatten()
    dRv = (-krep*Rv*Omrsv + krepr*Sv*ROmrsv + kdrd_ROmrsv)
    dSv = (krep*Rv*Omrsv - krepr*Sv*ROmrsv - kdrd_ROmrsv)
    dROm = (Om @ (krep*Rm) - ROm @ (krepr*Sm) - kdrd*ROm).flatten()
    
    duThv = (ktxnTh*Thtempv -krzTh*uThv)
    dThv = (krzTh*uThv -kth*Thv*Omrsv)
    
    dFv = (ktxnF*Ftempv - krsdF*Fv*GOmrsv + krevF_Omrsv*GFv - kssdF*Fv \
           -krsdF*Fv*AGObmrsv + krevF_Omrsv*AGFbv)
    dGFv = (krsdF*Fv*GOmrsv - krevF_Omrsv*GFv)
    
    duAGm = (ktxnAG*AGtempm - krzA*uAGm - kdsduAG*uAGm).flatten()
    dAGm = (krzA*uAGm - Omrsd @ (krsdA*AGm) - kdsdAG*AGm).flatten()
    dAGOam = (AGmap @ (Om @ (krsdA*AGm)) - Omrsd @ (krsdA*AGOam) + AGObmrsd @ (krevA*Om)).flatten()
    dAGObm = (Om @ krsdA_AGOamcsd -AGObm@krevA_Omcsd \
              -AGObm @ krsdFm + (krevF*Om) @ AGFbm).flatten()
    
    dAGFbv = krsdF*Fv*AGObmrsv - krevF_Omrsv*AGFbv
    
    duCGm = (ktxnCG*CGtempm - krzCG*uCGm - kdsduCG*uCGm).flatten()
    dCGm = (krzCG*uCGm  - (Omrsd @ (krsdCG*CGm).T).T + CGmap @ krevCG_CGObmrsd - Omrsd@(krsdCG*CGm) + krevCG_CGOamrsd @ CGmap - kdsdCG*CGm).flatten()
    dCGOam = (Om @ krsdCG_CGmcsd - krevCG*CGOam - (krsdCG*CGOam) @ CGannA).flatten()
    dCGObm = (Om @ krsdCG_CGmrsd - krevCG*CGObm - (krsdCG*CGObm) @ CGannB).flatten()
    
    
    
    return(np.concatenate([duGm,dGm,dOm,dGOm,dRv,dSv,dROm,duThv,dThv,dFv,dGFv,duAGm,dAGm,dAGOam,dAGObm,duCGm,dCGm,dCGOam,dCGObm,dAGFbv]))


class RSD_sim:
    def __init__(self,domains=5):
        
        self.N = domains
        # DNA template concentration vectors
        #self.IN_con = domains*[0] # Inputs
        
       
        self.Gtemp_con = 0*np.ones((domains,domains))  # single RNA gates
        self.Otemp_con = 0*np.ones((domains,domains)) # Outputs
        self.thresh_con = np.array(domains*[0])
        self.F_con = np.array(domains*[0])
        self.AG_con = 0*np.ones((domains,domains))
        self.CG_con = 0*np.ones((domains,domains))
        
        # Initial species concentration vectors
        self.uG_ic = 0*np.ones((domains,domains))
        self.GO_ic = 0*np.ones((domains,domains))
        self.G_ic = 0*np.ones((domains,domains)) # single RNA gates
        #self.rsdA_ic = domains*[0] # AND gates
        self.rep_ic = np.array(domains*[0]) # reporters
        self.RO_ic =  0*np.ones((domains,domains))
        self.S_ic = np.array(domains*[0])
        #self.f_ic = domains*[0] # fuel strands
        #self.wta_ic = domains*[0] # WTA strands
        #self.th_ic = domains*[0] # threshold gates
        self.out_ic = 0*np.ones((domains,domains)) # Outputs
        self.rep_ic_flag = np.array(domains*[0]) # saves if the reporter con was updated
        
        self.uTh_ic = np.array(domains*[0])
        self.thresh_ic = np.array(domains*[0])
        
        self.F_ic = np.array(domains*[0])
        self.GF_ic = np.array(domains*[0])
        
        self.uAG_ic = 0*np.ones((domains,domains))
        self.AG_ic = 0*np.ones((domains,domains))
        self.AGOa_ic = 0*np.ones((domains,domains))
        self.AGOb_ic = 0*np.ones((domains,domains))
        self.AGmap = 0*np.ones((domains,domains))
        
        self.uCG_ic = 0*np.ones((domains,domains))
        self.CG_ic = 0*np.ones((domains,domains))
        self.CGOa_ic = 0*np.ones((domains,domains))
        self.CGOb_ic = 0*np.ones((domains,domains))
        self.CGmap = 0*np.ones((domains,domains))  
        
        #self.AGFa_ic = 0*np.ones((domains,domains))
        self.AGFb_ic = np.array(domains*[0])
        
        self.ktxnO = 0.013*np.ones((domains,domains))
        self.ktxnG = 0.013*np.ones((domains,domains))
        self.ktxnTh = np.array(domains*[0.013])
        self.ktxnF = np.array(domains*[0.013])
        self.ktxnAG = 0.013*np.ones((domains,domains))
        self.ktxnCG = 0.013*np.ones((domains,domains))
        
        self.krz = (.25/60)*np.ones((domains,domains))
        self.krsd = (1e3/1e9)*np.ones((domains,domains))
        self.krev = (270/1e9)*np.ones((domains,domains))
         
        for x in range(self.N):
            self.krev[x,x] = 0
            
        self.krep = np.array(domains*[1e4 / 1e9])
        self.krepr = np.array(domains*[0])
        
        self.kth = np.array(domains*[1e5 / 1e9])
        self.krzTh = np.array(domains*[.00417])
        
        self.krsdF = np.array(domains*[1e3 / 1e9])
        #self.krevF = np.array(domains*[270 / 1e9])
        self.krevF = (1e3/1e9)*np.ones((domains,domains))
        
        self.krzA = (.25/60)*np.ones((domains,domains))
        self.krsdA = (1e3/1e9)*np.ones((domains,domains))
        self.krevA = (270/1e9)*np.ones((domains,domains))
         
        for x in range(self.N):
            self.krevA[x,x] = 0
            
        self.krzCG = .00417*np.ones((domains,domains))
        self.krsdCG = (1e5/1e9)*np.ones((domains,domains))
        self.krevCG = (1)*np.ones((domains,domains))
        
        self.kssdO = 0*np.ones((domains,domains))
        self.kssdF = np.array(domains*[0])
        
        self.kdsduG = 0*np.ones((domains,domains))
        self.kdsdG = 0*np.ones((domains,domains))
        self.kdsdGO = 0*np.ones((domains,domains))
        self.kdsduAG = 0*np.ones((domains,domains))
        self.kdsdAG = 0*np.ones((domains,domains))
        self.kdsduCG = 0*np.ones((domains,domains))
        self.kdsdCG = 0*np.ones((domains,domains))
        
        self.kdrd = 0*np.ones((domains,domains))
        
        self.AGcheck = False
        self.Fcheck = False    
        
        self.initialcheck = np.array((7*self.N+13*self.N**2)*[0])
        self.initialcheckIter = []
        
        self.plotCheck1 = False
        self.plotCheck2 = False
        
        
       
        
    # function for defining the DNA species in the model instance and the circuit connectivity
    def global_rate_constants(self,krz='False',krsd='False',krev='False',krep='False',krepr='False',kth='False',krzTh='False',krsdF='False',krevF='False',krevA='False',krsdA='False',krzA='False',krevCG='False',krsdCG='False',krzCG='False',ktxnO='False',ktxnG='False',ktxnTh='False',ktxnF='False',ktxnAG='False',ktxnCG='False',ktxn='False',kssdO='False',kssdF='False',kdsduG='False',kdsdG='False',kdsdGO='False',kdsduAG='False',kdsdAG='False',kdsduCG='False',kdsdCG='False',kdrd='False',kdeg='False',kssd='False',kdsd='False'):
        if krz != 'False':
            self.krz = krz*np.ones((self.N,self.N))
        if krsd != 'False':
            self.krsd = krsd*np.ones((self.N,self.N))
        if krev != 'False':
            self.krev = krev*np.ones((self.N,self.N))
            for x in range(self.N):
                self.krev[x,x] = 0
        if krep != 'False':
            self.krep = np.array(self.N*[krep])
        if krepr != 'False':
            self.krepr = np.array(self.N*[krepr])
        if kth != 'False':
            self.kth = np.array(self.N*[kth])
        if krzTh != 'False':
            self.krzTh = np.array(self.N*[krzTh])
        if krsdF != 'False':
            self.krsdF = np.array(self.N*[krsdF])
        if krevF != 'False':
            self.krevF = np.array(self.N*[krevF])
        if krsdA != 'False':
            self.krsdA = krsdA*np.ones((self.N,self.N))
        if krevA != 'False':
            self.krevA = krevA*np.ones((self.N,self.N))
            for x in range(self.N):
                self.krevA[x,x] = 0
        if krzA != 'False':
            self.krzA = krzA*np.ones((self.N,self.N))
        if krsdCG != 'False':
            self.krsdCG = krsdCG*np.ones((self.N,self.N))
        if krevCG != 'False':
            self.krevCG = krevCG*np.ones((self.N,self.N))
        if krzCG != 'False':
            self.krzCG = krzCG*np.ones((self.N,self.N))
        if ktxnO != 'False':
            self.ktxnO = ktxnO*np.ones((self.N,self.N))
        if ktxnG != 'False':
            self.ktxnG = ktxnG*np.ones((self.N,self.N))
        if ktxnTh != 'False':
            self.ktxnTh = np.array(self.N*[ktxnTh])
        if ktxnF != 'False':
            self.ktxnF = np.array(self.N*[ktxnF])
        if ktxnAG != 'False':
            self.ktxnAG = ktxnAG*np.ones((self.N,self.N))
        if ktxnCG != 'False':
            self.ktxnCG = ktxnCG*np.ones((self.N,self.N))
        if ktxn != 'False':
            self.ktxnO = ktxn*np.ones((self.N,self.N))
            self.ktxnG = ktxn*np.ones((self.N,self.N))
            self.ktxnTh = np.array(self.N*[ktxn])
            self.ktxnF = np.array(self.N*[ktxn])
            self.ktxnAG = ktxn*np.ones((self.N,self.N))
            self.ktxnCG = ktxn*np.ones((self.N,self.N))
        
        if kssdO != 'False':
            self.kssdO = kssdO*np.ones((self.N,self.N))
        if kssdF != 'False':
            self.kssdF = np.array(self.N*[kssdF])    
        if kdsduG != 'False':
            self.kdsduG = kdsduG*np.ones((self.N,self.N))
        if kdsdG != 'False':
            self.kdsdG = kdsdG*np.ones((self.N,self.N))
        if kdsdGO != 'False':
            self.kdsdGO = kdsdGO*np.ones((self.N,self.N))
        if kdsduAG != 'False':
            self.kdsduAG = kdsduAG*np.ones((self.N,self.N))
        if kdsdAG != 'False':
            self.kdsdAG = kdsdAG*np.ones((self.N,self.N))
        if kdsduCG != 'False':
            self.kdsduCG = kdsduCG*np.ones((self.N,self.N))
        if kdsdCG != 'False':
            self.kdsdCG = kdsdCG*np.ones((self.N,self.N))
        
        if kdrd != 'False':
            self.kdrd = kdrd*np.ones((self.N,self.N))
            
        if kdeg != 'False':
            self.kssdO = kdeg*np.ones((self.N,self.N))
            self.kssdF = np.array(self.N*[kdeg])
            self.kdsduG = kdeg*np.ones((self.N,self.N))
            self.kdsdG = kdeg*np.ones((self.N,self.N))
            self.kdsdGO = kdeg*np.ones((self.N,self.N))
            self.kdsduAG = kdeg*np.ones((self.N,self.N))
            self.kdsdAG = kdeg*np.ones((self.N,self.N))
            self.kdsdsuCG = kdeg*np.ones((self.N,self.N))
            self.kdsdCG = kdeg*np.ones((self.N,self.N))
            self.kdrd = kdeg*np.ones((self.N,self.N))
        
        if kssd != 'False':
            self.kssdO = kssd*np.ones((self.N,self.N))
            self.kssdF = np.array(self.N*[kssd])
        
        if kdsd != 'False':
            self.kdsduG = kdsd*np.ones((self.N,self.N))
            self.kdsdG = kdsd*np.ones((self.N,self.N))
            self.kdsdGO = kdsd*np.ones((self.N,self.N))
            self.kdsduAG = kdsd*np.ones((self.N,self.N))
            self.kdsdAG = kdsd*np.ones((self.N,self.N))
            self.kdsdsuCG = kdsd*np.ones((self.N,self.N))
            self.kdsdCG = kdsd*np.ones((self.N,self.N))
            
    
    
    def molecular_species(self,name,DNA_con=0,ic='False',krz='False',krsd='False',krev='False',krep='False',krepr='False',kth='False',krzTh='False',krsdF='False',krevF='False',krevA='False',krsdA='False',krzA='False',krevCG='False',krsdCG='False',krzCG='False',ktxnO='False',ktxnG='False',ktxnTh='False',ktxnF='False',ktxnAG='False',ktxnCG='False',kssdO='False',kssdF='False',kdsduG='False',kdsdG='False',kdsdGO='False',kdsduAG='False',kdsdAG='False',kdsduCG='False',kdsdCG='False',kdrd='False'):
        inp = re.compile("(i|in|inp|input)\{\w*\d+\w*\}")
        inps = inp.fullmatch(name.lower())
        rep = re.compile("(r|rep|reporter)\{\w*\d+\w*\}")
        reps = rep.fullmatch(name.lower())
        out = re.compile("(o|out|output)\{\w*\d+\w*\,\w*\d+\w*\}")
        outs = out.fullmatch(name.lower())
        gate = re.compile("(g|gate)\{\w*\d+\w*\,\w*\d+\w*\}")
        gates = gate.fullmatch(name.lower())
        
        uG = re.compile("ug\{\w*\d+\w*\,\w*\d+\w*\}")
        uGs = uG.fullmatch(name.lower())
        GI = re.compile("gi\{\w*\d+\w*\}")
        GIs = GI.fullmatch(name.lower())
        GO = re.compile("go\{\w*\d+\w*\,\w*\d+\w*\}")
        GOs = GO.fullmatch(name.lower())
        RO = re.compile("ro\{\w*\d+\w*\,\w*\d+\w*\}")
        ROs = RO.fullmatch(name.lower())
        S =  re.compile("s\{\w*\d+\w*\}")
        Ss = S.fullmatch(name.lower())
        
        uTh = re.compile("uth\{\w*\d+\w*\}")
        uThs = uTh.fullmatch(name.lower())
        thresh = re.compile("th\{\w*\d+\w*\}")
        threshs = thresh.fullmatch(name.lower())
        
        fuel = re.compile("f\{\w*\d+\w*\}")
        fuels = fuel.fullmatch(name.lower())
        GF = re.compile("f\{\w*\d+\w*\}")
        GFs = GF.fullmatch(name.lower())
        
        #AG = re.compile("ag\{\w*\d+(\.\d+)+\,\w*\d+\w*\}")
        uAG = re.compile("uag\{\w*\d+\w*\,\w*\d+\w*\}")
        uAGs = uAG.fullmatch(name.lower())
        AG = re.compile("(ag|g|gate)\{\w*\d+\w*(\.|\&)\w*\d+\w*\,\w*\d+\w*\}")
        AGs = AG.fullmatch(name.lower())
        AGOa = re.compile("agoa\{\w*\d+\w*\,\w*\d+\w*\}")
        AGOas = AGOa.fullmatch(name.lower())
        AGOb = re.compile("agob\{\w*\d+\w*\,\w*\d+\w*\}")
        AGObs = AGOb.fullmatch(name.lower())
        
        #AGFa  = re.compile("agfa\{\w*\d+\w*\,\w*\d+\w*\}")
        #AGFas = AGFa.fullmatch(name.lower())
        AGFb = re.compile("agfb\{\w*\d+\w*\}")
        AGFbs = AGFb.fullmatch(name.lower())
        
        uCG = re.compile("ucg\{\w*\d+\w*\,\w*\d+\w*\}")
        uCGs = uCG.fullmatch(name.lower())
        CG = re.compile("cg\{\w*\d+\w*\,\w*\d+\w*\}")
        CGs = CG.fullmatch(name.lower())
        CGOa = re.compile("cgoa\{\w*\d+\w*\,\w*\d+\w*\}")
        CGOas = CGOa.fullmatch(name.lower())
        CGOb = re.compile("cgob\{\w*\d+\w*\,\w*\d+\w*\}")
        CGObs = CGOb.fullmatch(name.lower())
        
        
        if inps:
            inpInd1f = re.compile("\d+")
            inpInd1 = int(inpInd1f.search(name.lower()).group())-1
                
            self.Otemp_con[inpInd1,inpInd1]=DNA_con
            if ic != 'False':
                self.out_ic[inpInd1,inpInd1]=ic
                self.initialcheck[2*self.N**2+(self.N*inpInd1) + inpInd1] += 1
                
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to an input')
            if ktxnO != 'False':
                self.ktxnO[inpInd1,inpInd1] = ktxnO
            if kssdO != 'False':
                self.kssdO[inpInd1,inpInd1] = kssdO
            
        elif reps:
            repInd1f = re.compile('\d+')
            repInd1 = int(repInd1f.search(name.lower()).group())-1
            
            if ic == 'False':
                self.rep_ic[repInd1]=DNA_con
            else:
                self.rep_ic[repInd1] = ic
                self.initialcheck[4*self.N**2 + repInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a reporter')
            
            if krep != 'False':
                self.krep[repInd1] = krep
            if krepr != 'False':
                self.krepr[repInd1] = krepr
        elif outs:
            outIndf = re.findall("\d+",name.lower())

            outInd1 = int(outIndf[0])-1
            outInd2 = int(outIndf[1])-1
            
            
            self.Otemp_con[outInd1,outInd2]=DNA_con
            if ic != 'False':
                self.out_ic[outInd1,outInd2]=ic
                self.initialcheck[2*self.N**2+(self.N*outInd1) + outInd2] += 1
            
            
            if krz != 'False' or krsd != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False'  or krzA != 'False' or krsdA != "False" or krsdCG != 'False' or krzCG != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to an output')
            if ktxnO != 'False':
                self.ktxnO[outInd1,outInd2] = ktxnO
            if krev != 'False':
                self.krev[outInd1,outInd2] = krev
            if krevF != 'False':
                self.krevF[outInd1,outInd2] = krevF
            if krevA != 'False':
                self.krevA[outInd1,outInd2] = krevA
            if krevCG != 'False':
                self.krevCG[outInd1,outInd2] = krevCG
            if kssdO != 'False':
                self.kssdO[outInd1,outInd2] = kssdO
            
            
            
        elif gates:
            gateIndf = re.findall("\d+",name.lower())
            gateInd1 = int(gateIndf[0])-1
            gateInd2 = int(gateIndf[1])-1
            
            self.Gtemp_con[gateInd1,gateInd2]=DNA_con
            if ic != 'False':
                self.G_ic[gateInd1,gateInd2]=ic
                self.initialcheck[self.N**2+(self.N*gateInd1) + gateInd2] += 1
            
            if krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a gate')
            if krev != 'False':
                self.krev[gateInd1,gateInd2] = krev   
            if krsd != 'False':
                self.krsd[gateInd1,gateInd2] = krsd
            if krz != 'False':
                self.krz[gateInd1,gateInd2] = krz
            if ktxnG != 'False':
                self.ktxnG[gateInd1,gateInd2] = ktxnG
            if kdsdG != 'False':
                self.kdsdG[gateInd1,gateInd2] = kdsdG
            if krsdA != 'False':
                self.krsdA[gateInd1,gateInd2] = krsdA
    
        
        elif uGs:
            uGIndf = re.findall("\d+",name.lower())
            uGInd1 = int(uGIndf[0])-1
            uGInd2 = int(uGIndf[1])-1
            
            if ic != 'False':
                self.uG_ic[uGInd1,uGInd2]=ic
                self.initialcheck[0+(self.N*uGInd1) + uGInd2] += 1
            
            if krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to an uncleaved gate')
            if krz != 'False':
                self.krz[uGInd1,uGInd2] = krz
            if ktxnG != 'False':
                self.ktxnG[uGInd1,uGInd2] = ktxnG
            if kdsduG != 'False':
                self.kdsduG[uGInd1,uGInd2] = kdsduG
            
            
        elif GIs:
            GIInd1f = re.compile("\d+")
            GIInd1 = int(GIInd1f.search(name.lower()).group())-1
            
            if ic != 'False':
                self.GO_ic[GIInd1,GIInd1]=ic
                self.initialcheck[3*self.N**2+(self.N*GIInd1) + GIInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a GI species')
            if kdsdGO != 'False':
                self.kdsdGO[GIInd1,GIInd1] = kdsdGO
            
            
        elif GOs:
            GOIndf = re.findall("\d+",name.lower())
            GOInd1 = int(GOIndf[0])-1
            GOInd2 = int(GOIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to GO species')
            if kdsdGO != 'False':
                self.kdsdGO[GOInd1,GOInd2] = kdsdGO
            
            if ic != 'False':
                self.GO_ic[GOInd1,GOInd2]=ic
                self.initialcheck[3*self.N**2+(self.N*GOInd1) + GOInd2] += 1
            
        elif ROs:
           ROIndf = re.findall("\d+",name.lower())
           ROInd1 = int(ROIndf[0])-1
           ROInd2 = int(ROIndf[1])-1
           
           if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False':
               print('This rate constant should not be changed with respect to an RO species')
           if kdrd != 'False':
               self.kdrd[ROInd1,ROInd2] = kdrd
          
           if ic != 'False':
               self.RO_ic[ROInd1,ROInd2]=ic
               self.initialcheck[2*self.N+4*self.N**2+(self.N*ROInd1) + ROInd2] += 1
            
        elif Ss:
            SInd1f = re.compile("\d+")
            SInd1 = int(SInd1f.search(name.lower()).group())-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            
            if ic != 'False':
                self.S_ic[SInd1]=ic
                self.initialcheck[self.N+4*self.N**2 + SInd1] += 1
        
        elif uThs:
            uThInd1f = re.compile("\d+")
            uThInd1 = int(uThInd1f.search(name.lower()).group())-1
            
            if ic != 'False':
                self.uTh_ic[uThInd1] = ic
                self.initialcheck[2*self.N+5*self.N**2 + uThInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krzTh != 'False':
                self.krzTh[uThInd1] = krzTh
            if ktxnTh != 'False':
                self.ktxnTh[uThInd1] = ktxnTh
        
        elif threshs:
            threshInd1f = re.compile("\d+")
            threshInd1 = int(threshInd1f.search(name.lower()).group())-1
            
            self.thresh_con[threshInd1] = DNA_con
            if ic != 'False':
                self.thresh_ic[threshInd1] = ic
                self.initialcheck[3*self.N+5*self.N**2 + threshInd1] += 1
            
            if krep != 'False' or krepr != 'False' or krev != 'False' or krsd != 'False' or krz != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a gate')
            if kth != 'False':
                self.kth[threshInd1] = kth
            if krzTh != 'False':
                self.krzTh[threshInd1] = krzTh
            if ktxnTh != 'False':
                self.ktxnTh[threshInd1] = ktxnTh
        
        elif fuels:
            fuelInd1f = re.compile("\d+")
            self.fuelInd1 = int(fuelInd1f.search(name.lower()).group())-1
            
            self.F_con[self.fuelInd1] = DNA_con
            if ic != 'False':
                self.F_ic[self.fuelInd1] = ic
                self.initialcheck[4*self.N+5*self.N**2 + self.fuelInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krsdF != 'False':
                self.krsdF[self.fuelInd1] = krsdF
            if krevF != 'False':
                self.krevF[self.fuelInd1] = krevF
            if ktxnF != 'False':
                self.ktxnF[self.fuelInd1] = ktxnF
            if kssdF != 'False':
                self.kssdF[self.fuelInd1] = kssdF
            
            self.Fcheck = True
            
            
        elif GFs:
            GFInd1f = re.compile("\d+")
            GFInd1 = int(GFInd1f.search(name.lower()).group())-1
            
            self.GF_con[GFInd1] = DNA_con
            if ic != 'False':
                self.GF_ic[GFInd1] = ic
                self.initialcheck[5*self.N+5*self.N**2 + self.GFInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if ktxnF != 'False':
                self.ktxnF[GFInd1] = ktxnF
                
        elif AGs:
            AGIndf = re.findall("\d+",name.lower())
            self.AGInd1 = int(AGIndf[0])-1
            AGInd2 = int(AGIndf[-1])-1
            
            AGOInd1 = int(AGIndf[0])-1
            AGOInd2 = int(AGIndf[1])-1
            
            self.AG_con[self.AGInd1,AGInd2] = DNA_con
            if ic != 'False':
                self.AG_ic[self.AGInd1,AGInd2] = ic
                self.initialcheck[6*self.N+6*self.N**2:+(self.N*self.AGInd1) + AGInd2] += 1
            
            self.AGmap[AGOInd1,AGOInd2] = 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krsdA != 'False':
                self.krsdA[self.AGInd1,AGInd2] = krsdA
            if krevA != 'False':
                self.krevA[AGOInd2,AGInd2] = krevA
            if krzA != 'False':
                self.krzA[self.AGInd1,AGInd2] = krzA
            if ktxnAG != 'False':
                self.ktxnAG[self.AGInd1,AGInd2] = ktxnAG
            if kdsdAG != 'False':
                self.kdsdAG[self.AGInd1,AGInd2] = kdsdAG
                
            self.AGcheck = True
                
        elif uAGs:
            uAGIndf = re.findall("\d+",name.lower())
            uAGInd1 = int(uAGIndf[0])-1
            uAGInd2 = int(uAGIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krzA != 'False':
                self.krzA[uAGInd1,uAGInd2] = krzA
            if ktxnAG != 'False':
                self.ktxnAG[uAGInd1,uAGInd2] = ktxnAG
            if kdsduAG != 'False':
                self.kdsduAG[uAGInd1,uAGInd2] = kdsduAG
            
            if ic != 'False':
                self.uAG_ic[uAGInd1,uAGInd2] = ic
                self.initialcheck[6*self.N+5*self.N**2:+(self.N*uAGInd1) + uAGInd2] += 1
                
        elif AGOas:
            AGOaIndf = re.findall("\d+",name.lower())
            AGOaInd1 = int(AGOaIndf[0])-1
            AGOaInd2 = int(AGOaIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krsdA != 'False':
                self.krsdA[AGOaInd1,AGOaInd2] = krsdA
            
            if ic != 'False':
                self.AGOa_ic[AGOaInd1,AGOaInd2] = ic
                self.initialcheck[6*self.N+7*self.N**2:+(self.N*AGOaInd1) + AGOaInd2] += 1
                
        elif AGObs:
            AGObIndf = re.findall("\d+",name.lower())
            AGObInd1 = int(AGObIndf[0])-1
            AGObInd2 = int(AGObIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            
            if ic != 'False':
                self.AGOb_ic[AGObInd1,AGObInd2] = ic
                self.initialcheck[6*self.N+8*self.N**2:+(self.N*AGObInd1) + AGObInd2] += 1
        
        
        elif AGFbs:
            AGFbIndf = re.findall("\d+",name.lower())
            AGFbInd1 = int(AGFbIndf[0])-1
            
            if ic != 'False':
                self.AGFb_ic[AGFbInd1] = ic
                self.initialcheck[6*self.N+13*self.N**2 + self.AGFbInd1] += 1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
                   
        elif CGs:
            CGIndf = re.findall("\d+",name.lower())
            
            CGInd1 = int(CGIndf[0])-1
            CGInd2 = int(CGIndf[1])-1
            
            
            self.CG_con[CGInd1,CGInd2] = DNA_con
            if ic != 'False':
                self.CG_ic[CGInd1,CGInd2] = ic
                self.initialcheck[6*self.N+10*self.N**2:+(self.N*CGInd1) + CGInd2] += 1
            
            self.CGmap[CGInd1,CGInd2] = 1
            
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krsdCG != 'False':
                self.krsdCG[CGInd1,CGInd2] = krsdCG
            if krevCG != 'False':
                self.krevCG[CGInd1,CGInd2] = krevCG
            if krzCG != 'False':
                self.krzCG[CGInd1,CGInd2] = krzCG
            if ktxnCG != 'False':
                self.ktxnCG[CGInd1,CGInd2] = ktxnCG
            if kdsdCG != 'False':
                self.kdsdCG[CGInd1,CGInd2] = kdsdCG
                
        elif uCGs:
            uCGIndf = re.findall("\d+",name.lower())
            uCGInd1 = int(uCGIndf[0])-1
            uCGInd2 = int(uCGIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krzCG != 'False':
                self.krzCG[uCGInd1,uCGInd2] = krzCG
            if ktxnCG != 'False':
                self.ktxnCG[uCGInd1,uCGInd2] = ktxnCG
            if kdsduCG != 'False':
                self.kdsduCG[uCGInd1,uCGInd2] = kdsduCG
            
            if ic != 'False':
                self.uCG_ic[uCGInd1,uCGInd2] = ic
                self.initialcheck[6*self.N+9*self.N**2:+(self.N*uCGInd1) + uCGInd2] += 1
                
        elif CGOas:
            CGOaIndf = re.findall("\d+",name.lower())
            CGOaInd1 = int(CGOaIndf[0])-1
            CGOaInd2 = int(CGOaIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            if krsdCG != 'False':
                self.krsdCG[CGOaInd1,CGOaInd2] = krsdCG
            
            if ic != 'False':
                self.CGOa_ic[CGOaInd1,CGOaInd2] = ic
                self.initialcheck[6*self.N+11*self.N**2:+(self.N*CGOaInd1) + CGOaInd2] += 1
            
        elif CGObs:
            CGObIndf = re.findall("\d+",name.lower())
            CGObInd1 = int(CGObIndf[0])-1
            CGObInd2 = int(CGObIndf[1])-1
            
            if krz != 'False' or krsd != 'False' or krev != 'False' or krep != 'False' or krepr != 'False' or kth != 'False' or krzTh != 'False' or krsdF != 'False' or krevF != 'False' or krevA != 'False' or krsdA != 'False' or krzA != 'False' or krevCG != 'False' or krsdCG != 'False' or krzCG != 'False' or ktxnO != 'False' or ktxnG != 'False' or ktxnTh != 'False' or ktxnF != 'False' or ktxnAG != 'False' or ktxnCG != 'False' or kssdO != 'False' or kssdF != 'False' or kdsduG != "False" or kdsdG != 'False' or kdsdGO != 'False' or kdsduAG != 'False' or kdsdAG != 'False' or kdsduCG != 'False' or kdsdCG != 'False' or kdrd != 'False':
                print('This rate constant should not be changed with respect to a S species')
            
            if ic != 'False':
                self.CGOb_ic[CGObInd1,CGObInd2] = ic
                self.initialcheck[6*self.N+12*self.N**2:+(self.N*CGObInd1) + CGObInd2] += 1
                    
            
        else:
            print("Please check README for correct nomenclature")
        
        if (self.AGcheck and self.Fcheck):
            if (self.fuelInd1 == self.AGInd1):
                print('Warning: Fuel reactions do not function properly with first index of AND Gate')
        
                
    # function for simulating the model instance
    def simulate(self,t_vec,leak=0.03,leakA=0.06,smethod='False',iteration=1):
        
        
        if iteration == 1:
        
            self.initials = np.concatenate([self.uG_ic.flatten(),self.G_ic.flatten(),self.out_ic.flatten(),self.GO_ic.flatten(),self.rep_ic.flatten(),self.S_ic.flatten(),self.RO_ic.flatten(),self.uTh_ic.flatten(),self.thresh_ic.flatten(),self.F_ic.flatten(),self.GF_ic.flatten(),self.uAG_ic.flatten(),self.AG_ic.flatten(),self.AGOa_ic.flatten(),self.AGOb_ic.flatten(),self.uCG_ic.flatten(),self.CG_ic.flatten(),self.CGOa_ic.flatten(),self.CGOb_ic.flatten(),self.AGFb_ic.flatten()])
            self.initialcheckIter.append(np.asarray(list(self.initialcheck)))
            
            if smethod != 'False':
                self.solve = spi.solve_ivp(lambda t,y: rate_eqs(t,y,self.ktxnO,self.ktxnG,self.ktxnTh,self.ktxnF,self.ktxnAG,self.ktxnCG,self.krz,self.krsd,self.krev,self.krep,self.krepr,self.kth,self.krzTh,self.krsdF,self.krevF,self.krevA,self.krsdA,self.krzA,self.krevCG,self.krsdCG,self.krzCG,leak,leakA,self.Otemp_con,self.Gtemp_con,self.thresh_con,self.F_con,self.AG_con,self.AGmap.T,self.CG_con,self.CGmap,self.kssdO,self.kssdF,self.kdsduG,self.kdsdG,self.kdsdGO,self.kdsduAG,self.kdsdAG,self.kdsduCG,self.kdsdCG,self.kdrd,self.N),[t_vec[0],t_vec[-1]],self.initials,t_eval=t_vec,method=smethod)
            else:
                self.solve = spi.solve_ivp(lambda t,y: rate_eqs(t,y,self.ktxnO,self.ktxnG,self.ktxnTh,self.ktxnF,self.ktxnAG,self.ktxnCG,self.krz,self.krsd,self.krev,self.krep,self.krepr,self.kth,self.krzTh,self.krsdF,self.krevF,self.krevA,self.krsdA,self.krzA,self.krevCG,self.krsdCG,self.krzCG,leak,leakA,self.Otemp_con,self.Gtemp_con,self.thresh_con,self.F_con,self.AG_con,self.AGmap.T,self.CG_con,self.CGmap,self.kssdO,self.kssdF,self.kdsduG,self.kdsdG,self.kdsdGO,self.kdsduAG,self.kdsdAG,self.kdsduCG,self.kdsdCG,self.kdrd,self.N),[t_vec[0],t_vec[-1]],self.initials,t_eval=t_vec,method='LSODA')
            
            
            self.t = self.solve.t
            self.uG_simcon = self.solve.y[0:self.N**2]
            self.G_simcon = self.solve.y[self.N**2:2*self.N**2]
            self.O_simcon = self.solve.y[2*self.N**2:3*self.N**2]
            self.GO_simcon = self.solve.y[3*self.N**2:4*self.N**2]
            self.R_simcon = self.solve.y[4*self.N**2:self.N+4*self.N**2]
            self.S_simcon = self.solve.y[self.N+4*self.N**2:2*self.N+4*self.N**2]
            self.RO_simcon = self.solve.y[2*self.N+4*self.N**2:2*self.N+5*self.N**2]
            self.uTh_simcon = self.solve.y[2*self.N+5*self.N**2:3*self.N+5*self.N**2]
            self.Th_simcon = self.solve.y[3*self.N+5*self.N**2:4*self.N+5*self.N**2]
            self.F_simcon = self.solve.y[4*self.N+5*self.N**2:5*self.N+5*self.N**2]
            self.GF_simcon = self.solve.y[5*self.N+5*self.N**2:6*self.N+5*self.N**2]
            self.uAG_simcon = self.solve.y[6*self.N+5*self.N**2:6*self.N+6*self.N**2]
            self.AG_simcon = self.solve.y[6*self.N+6*self.N**2:6*self.N+7*self.N**2]
            self.AGOa_simcon = self.solve.y[6*self.N+7*self.N**2:6*self.N+8*self.N**2]
            self.AGOb_simcon = self.solve.y[6*self.N+8*self.N**2:6*self.N+9*self.N**2]
            self.uCG_simcon = self.solve.y[6*self.N+9*self.N**2:6*self.N+10*self.N**2]
            self.CG_simcon = self.solve.y[6*self.N+10*self.N**2:6*self.N+11*self.N**2]
            self.CGOa_simcon = self.solve.y[6*self.N+11*self.N**2:6*self.N+12*self.N**2]
            self.CGOb_simcon = self.solve.y[6*self.N+12*self.N**2:6*self.N+13*self.N**2]
            self.AGFb_simcon = self.solve.y[6*self.N+13*self.N**2:7*self.N+13*self.N**2]
        
        if iteration > 1:
                     
            self.initialsNEW = np.array(self.solve.y[:,-1])
            self.initialsOLD = np.concatenate([self.uG_ic.flatten(),self.G_ic.flatten(),self.out_ic.flatten(),self.GO_ic.flatten(),self.rep_ic.flatten(),self.S_ic.flatten(),self.RO_ic.flatten(),self.uTh_ic.flatten(),self.thresh_ic.flatten(),self.F_ic.flatten(),self.GF_ic.flatten(),self.uAG_ic.flatten(),self.AG_ic.flatten(),self.AGOa_ic.flatten(),self.AGOb_ic.flatten(),self.uCG_ic.flatten(),self.CG_ic.flatten(),self.CGOa_ic.flatten(),self.CGOb_ic.flatten(),self.AGFb_ic.flatten()])
            self.initialcheckIter.append(np.asarray(list(self.initialcheck)))
            
            
            
            for i in range(len(self.initialsOLD)):
                if self.initialcheckIter[iteration-1][i] != self.initialcheckIter[iteration-2][i]:
                    self.initialsNEW[i] = self.initialsOLD[i]
            
            
            
            if smethod != 'False':
                self.solveNEW = spi.solve_ivp(lambda t,y: rate_eqs(t,y,self.ktxnO,self.ktxnG,self.ktxnTh,self.ktxnF,self.ktxnAG,self.ktxnCG,self.krz,self.krsd,self.krev,self.krep,self.krepr,self.kth,self.krzTh,self.krsdF,self.krevF,self.krevA,self.krsdA,self.krzA,self.krevCG,self.krsdCG,self.krzCG,leak,leakA,self.Otemp_con,self.Gtemp_con,self.thresh_con,self.F_con,self.AG_con,self.AGmap.T,self.CG_con,self.CGmap,self.kssdO,self.kssdF,self.kdsduG,self.kdsdG,self.kdsdGO,self.kdsduAG,self.kdsdAG,self.kdsduCG,self.kdsdCG,self.kdrd,self.N),[t_vec[0],t_vec[-1]],self.initialsNEW,t_eval=t_vec,method=smethod)
            else:
                self.solveNEW = spi.solve_ivp(lambda t,y: rate_eqs(t,y,self.ktxnO,self.ktxnG,self.ktxnTh,self.ktxnF,self.ktxnAG,self.ktxnCG,self.krz,self.krsd,self.krev,self.krep,self.krepr,self.kth,self.krzTh,self.krsdF,self.krevF,self.krevA,self.krsdA,self.krzA,self.krevCG,self.krsdCG,self.krzCG,leak,leakA,self.Otemp_con,self.Gtemp_con,self.thresh_con,self.F_con,self.AG_con,self.AGmap.T,self.CG_con,self.CGmap,self.kssdO,self.kssdF,self.kdsduG,self.kdsdG,self.kdsdGO,self.kdsduAG,self.kdsdAG,self.kdsduCG,self.kdsdCG,self.kdrd,self.N),[t_vec[0],t_vec[-1]],self.initialsNEW,t_eval=t_vec,method='LSODA')
            
            self.solve.t = np.concatenate([self.solve.t,self.solveNEW.t])
            self.solve.y = np.concatenate([self.solve.y,self.solveNEW.y],axis=1)
            
            #self.solve.t = self.solveNEW.t
            #self.solve.y = self.solveNEW.y
            
            self.t = self.solve.t
            self.uG_simcon = self.solve.y[0:self.N**2]
            self.G_simcon = self.solve.y[self.N**2:2*self.N**2]
            self.O_simcon = self.solve.y[2*self.N**2:3*self.N**2]
            self.GO_simcon = self.solve.y[3*self.N**2:4*self.N**2]
            self.R_simcon = self.solve.y[4*self.N**2:self.N+4*self.N**2]
            self.S_simcon = self.solve.y[self.N+4*self.N**2:2*self.N+4*self.N**2]
            self.RO_simcon = self.solve.y[2*self.N+4*self.N**2:2*self.N+5*self.N**2]
            self.uTh_simcon = self.solve.y[2*self.N+5*self.N**2:3*self.N+5*self.N**2]
            self.Th_simcon = self.solve.y[3*self.N+5*self.N**2:4*self.N+5*self.N**2]
            self.F_simcon = self.solve.y[4*self.N+5*self.N**2:5*self.N+5*self.N**2]
            self.GF_simcon = self.solve.y[5*self.N+5*self.N**2:6*self.N+5*self.N**2]
            self.uAG_simcon = self.solve.y[6*self.N+5*self.N**2:6*self.N+6*self.N**2]
            self.AG_simcon = self.solve.y[6*self.N+6*self.N**2:6*self.N+7*self.N**2]
            self.AGOa_simcon = self.solve.y[6*self.N+7*self.N**2:6*self.N+8*self.N**2]
            self.AGOb_simcon = self.solve.y[6*self.N+8*self.N**2:6*self.N+9*self.N**2]
            self.uCG_simcon = self.solve.y[6*self.N+9*self.N**2:6*self.N+10*self.N**2]
            self.CG_simcon = self.solve.y[6*self.N+10*self.N**2:6*self.N+11*self.N**2]
            self.CGOa_simcon = self.solve.y[6*self.N+11*self.N**2:6*self.N+12*self.N**2]
            self.CGOb_simcon = self.solve.y[6*self.N+12*self.N**2:6*self.N+13*self.N**2]
            self.AGFb_simcon = self.solve.y[6*self.N+13*self.N**2:7*self.N+13*self.N**2]
        
    def output_concentration(self,name):
        inp = re.compile("(i|in|inp|input)\{\w*\d+\w*\}")
        inps = inp.fullmatch(name.lower())
        rep = re.compile("(r|rep|reporter)\{\w*\d+\w*\}")
        reps = rep.fullmatch(name.lower())
        out = re.compile("(o|out|output)\{\w*\d+\w*\,\w*\d+\w*\}")
        outs = out.fullmatch(name.lower())
        gate = re.compile("(g|gate)\{\w*\d+\w*\,\w*\d+\w*\}")
        gates = gate.fullmatch(name.lower())
        
        uG = re.compile("ug\{\w*\d+\w*\,\w*\d+\w*\}")
        uGs = uG.fullmatch(name.lower())
        GI = re.compile("gi\{\w*\d+\w*\}")
        GIs = GI.fullmatch(name.lower())
        GO = re.compile("go\{\w*\d+\w*\,\w*\d+\w*\}")
        GOs = GO.fullmatch(name.lower())
        RO = re.compile("ro\{\w*\d+\w*\,\w*\d+\w*\}")
        ROs = RO.fullmatch(name.lower())
        S =  re.compile("s\{\w*\d+\w*\}")
        Ss = S.fullmatch(name.lower())
        
        uTh = re.compile("uth\{\w*\d+\w*\}")
        uThs = uTh.fullmatch(name.lower())
        thresh = re.compile("th\{\w*\d+\w*\}")
        threshs = thresh.fullmatch(name.lower())
        
        fuel = re.compile("f\{\w*\d+\w*\}")
        fuels = fuel.fullmatch(name.lower())
        GF = re.compile("f\{\w*\d+\w*\}")
        GFs = GF.fullmatch(name.lower())
        
        #AG = re.compile("ag\{\w*\d+(\.\d+)+\,\w*\d+\w*\}")
        uAG = re.compile("uag\{\w*\d+\w*\,\w*\d+\w*\}")
        uAGs = uAG.fullmatch(name.lower())
        AG = re.compile("(ag|g|gate)\{\w*\d+\w*(\.|\&)\w*\d+\w*\,\w*\d+\w*\}")
        AGs = AG.fullmatch(name.lower())
        AGOa = re.compile("agoa\{\w*\d+\w*\,\w*\d+\w*\}")
        AGOas = AGOa.fullmatch(name.lower())
        AGOb = re.compile("agob\{\w*\d+\w*\,\w*\d+\w*\}")
        AGObs = AGOb.fullmatch(name.lower())
        
        #AGFa  = re.compile("agfa\{\w*\d+\w*\,\w*\d+\w*\}")
        #AGFas = AGFa.fullmatch(name.lower())
        AGFb = re.compile("agfb\{\w*\d+\w*\}")
        AGFbs = AGFb.fullmatch(name.lower())
        
        uCG = re.compile("ucg\{\w*\d+\w*\,\w*\d+\w*\}")
        uCGs = uCG.fullmatch(name.lower())
        CG = re.compile("cg\{\w*\d+\w*\,\w*\d+\w*\}")
        CGs = CG.fullmatch(name.lower())
        CGOa = re.compile("cgoa\{\w*\d+\w*\,\w*\d+\w*\}")
        CGOas = CGOa.fullmatch(name.lower())
        CGOb = re.compile("cgob\{\w*\d+\w*\,\w*\d+\w*\}")
        CGObs = CGOb.fullmatch(name.lower())
        
        if inps:
            inpInd1f = re.compile("\d+")
            inpInd1 = int(inpInd1f.search(name.lower()).group())-1
            
            return(self.O_simcon[self.N*inpInd1 + inpInd1])
                
            
        elif reps:
            repInd1f = re.compile('\d+')
            repInd1 = int(repInd1f.search(name.lower()).group())-1
            
            return(self.R_simcon[repInd1])
        elif outs:
            outIndf = re.findall("\d+",name.lower())

            outInd1 = int(outIndf[0])-1
            outInd2 = int(outIndf[1])-1
            
            return(self.O_simcon[self.N*outInd1 + outInd2])
        
            
        elif gates:
            gateIndf = re.findall("\d+",name.lower())

            gateInd1 = int(gateIndf[0])-1
            gateInd2 = int(gateIndf[1])-1
            
            return(self.G_simcon[self.N*gateInd1 + gateInd2])
    
        
        elif uGs:
            uGIndf = re.findall("\d+",name.lower())
            uGInd1 = int(uGIndf[0])-1
            uGInd2 = int(uGIndf[1])-1
            
            return(self.uG_simcon[self.N*uGInd1 + uGInd2])        
            
        elif GIs:
            GIInd1f = re.compile("\d+")
            GIInd1 = int(GIInd1f.search(name.lower()).group())-1
            
            return(self.GO_simcon[self.N*GIInd1 + GIInd1])
        elif GOs:
            GOIndf = re.findall("\d+",name.lower())
            GOInd1 = int(GOIndf[0])-1
            GOInd2 = int(GOIndf[1])-1
            
            return(self.GO_simcon[self.N*GOInd1 + GOInd2])
            
        elif ROs:
            ROIndf = re.findall("\d+",name.lower())
            ROInd1 = int(ROIndf[0])-1
            ROInd2 = int(ROIndf[1])-1
            
            
            return(self.RO_simcon[self.N*ROInd1 + ROInd2])
        elif Ss:
            SInd1f = re.compile("\d+")
            SInd1 = int(SInd1f.search(name.lower()).group())-1
                
            return(self.S_simcon[SInd1])
        
        elif uThs:
            uThInd1f = re.compile("\d+")
            uThInd1 = int(uThInd1f.search(name.lower()).group())-1
                     
            return(self.uTh_simcon[uThInd1])
        
        elif threshs:
            threshInd1f = re.compile("\d+")
            threshInd1 = int(threshInd1f.search(name.lower()).group())-1
               
            return(self.Th_simcon[threshInd1])
        
        elif fuels:
            fuelInd1f = re.compile("\d+")
            fuelInd1 = int(fuelInd1f.search(name.lower()).group())-1
            
            return(self.F_simcon[fuelInd1])
                
        elif GFs:
            GFInd1f = re.compile("\d+")
            GFInd1 = int(GFInd1f.search(name.lower()).group())-1
            
            return(self.GF_simcon[GFInd1]) 
        elif AGs:
            AGIndf = re.findall("\d+",name.lower())
            AGInd1 = int(AGIndf[0])-1
            AGInd2 = int(AGIndf[-1])-1
            
            
            return(self.AG_simcon[self.N*AGInd1 + AGInd2])
        
        elif uAGs:
            uAGIndf = re.findall("\d+",name.lower())
            uAGInd1 = int(uAGIndf[0])-1
            uAGInd2 = int(uAGIndf[1])-1
            
            return(self.uAG_simcon[self.N*uAGInd1 + uAGInd2])
        
        elif AGOas:
            AGOaIndf = re.findall("\d+",name.lower())
            AGOaInd1 = int(AGOaIndf[0])-1
            AGOaInd2 = int(AGOaIndf[1])-1
            
            return(self.AGOa_simcon[self.N*AGOaInd1 + AGOaInd2])
        
        elif AGObs:
            AGObIndf = re.findall("\d+",name.lower())
            AGObInd1 = int(AGObIndf[0])-1
            AGObInd2 = int(AGObIndf[1])-1
            
            return(self.AGOb_simcon[self.N*AGObInd1 + AGObInd2])
        
        elif AGFbs:
            AGFbIndf = re.findall("\d+",name.lower())
            AGFbInd1 = int(AGFbIndf[0])-1
            
            
            return(self.AGFb_simcon[AGFbInd1])
        
        elif CGs:
            CGIndf = re.findall("\d+",name.lower())
            
            CGInd1 = int(CGIndf[0])-1
            CGInd2 = int(CGIndf[1])-1
           
            return(self.CG_simcon[self.N*CGInd1 + CGInd2])
            
            
                
        elif uCGs:
            uCGIndf = re.findall("\d+",name.lower())
            uCGInd1 = int(uCGIndf[0])-1
            uCGInd2 = int(uCGIndf[1])-1
            
            return(self.uCG_simcon[self.N*uCGInd1 + uCGInd2])
            
            
            
        elif CGOas:
            CGOaIndf = re.findall("\d+",name.lower())
            CGOaInd1 = int(CGOaIndf[0])-1
            CGOaInd2 = int(CGOaIndf[1])-1
            
            return(self.CGOa_simcon[self.N*CGOaInd1 + CGOaInd2])
            
            
            
        elif CGObs:
            CGObIndf = re.findall("\d+",name.lower())
            CGObInd1 = int(CGObIndf[0])-1
            CGObInd2 = int(CGObIndf[1])-1
            
            return(self.CGOb_simcon[self.N*CGObInd1 + CGObInd2])
            
            
        
        else:
            print("Please check README for correct nomenclature")
            
            
    def transcription_calibration(self,simTime,data,ktxn='False'):
        
        
        
        if ktxn != 'False':
            
            self.plotCheck1 = True
            
            if self.plotCheck2 == True:
                plt.figure()
            else:
                plt.figure(1)
            
            
            
            import Simulatorv2021 as RSDs #import version 2.0.2.1 of simulator

            fs = 12 #fontsize

            # Modeling
            if type(ktxn) != list:
                k_txn = [ktxn] #transcription rates
            else:
                k_txn = ktxn

            REP_con = 500 #reporter concentration

            model = RSDs.RSD_sim(5) # define the model instance and # of domains

            # initialize species involved in the reaction
            model.molecular_species('O{1,2}',DNA_con=25) 
            model.molecular_species('REP{2}',DNA_con=REP_con)
            
            
            for n in range(len(k_txn)):           
            
                model.global_rate_constants(ktxn=k_txn[n]) #globally changes transcription rates
                
                # run simulaton (input is simulaton time)
                model.simulate(simTime,smethod='LSODA')
                
                # pull out the species from the model solution to plot
                S2 = model.output_concentration('S{2}')
                    
            
                if len(k_txn) <= 4:
                    plt.subplot(2,4,1+n)
                    if (1+n) == 1:
                        plt.ylabel('Reacted reporter (%)',fontsize=fs)
                    ax1 = plt.gca()
                    ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
                    ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
                else:
                    if (1+n) <= (len(k_txn) / 2):
                        plt.subplot(3,int(math.ceil(len(k_txn)/2)),1+n)
                    else:
                        plt.subplot(3,int(math.ceil(len(k_txn)/2)),int(math.ceil((n+1)-len(k_txn)/2))+2*int(math.ceil(len(k_txn)/2)))   
                    ax1 = plt.gca()
                    ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
                    plt.tick_params(labelleft=False)
                    if (1+n) == 1 or (1+int(math.ceil((n+1)-len(k_txn)/2))+2*int(math.ceil(len(k_txn)/2))) % (int(math.ceil(len(k_txn)/2)) + 1) == 0:
                        plt.ylabel('Reacted reporter (%)',fontsize=fs)
                        plt.tick_params(labelleft=True, left=True)
                        ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
                    
   
                plt.plot(model.t/60,(data/REP_con)*100,color='aqua',linewidth=2,linestyle='-')
                plt.plot(model.t/60,(S2/REP_con)*100,color='orange',linewidth=2,linestyle='--')
                plt.xticks(fontsize=fs)
                plt.yticks(fontsize=fs)
                plt.ylim(-10,110)
                plt.xlim(0,120)
                ax1 = plt.gca()
                ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
                ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
                plt.xlabel('Time (min)',fontsize=fs)    
                plt.title('k='+str(k_txn[n]))
                plt.legend(['Raw Data','Calibrated'],frameon=False)
                
            
            
        
        else:
            
            self.plotCheck2 = True
            
            if self.plotCheck1 == True:
                plt.figure()
            else:
                plt.figure(1)
            
            
            
            import Simulatorv2021 as RSDs #import version 2.0.2.1 of simulator

            fs = 12 #fontsize

            # Modeling

            k_txn = [0.005,0.0075,0.01,0.0125,.015,.02] #transcription rates

            REP_con = 500 #reporter concentration

            model = RSDs.RSD_sim(5) # define the model instance and # of domains

            # initialize species involved in the reaction
            model.molecular_species('O{1,2}',DNA_con=25) 
            model.molecular_species('REP{2}',DNA_con=REP_con)
             
            
            for n in range(len(k_txn)):
                
                model.global_rate_constants(ktxn=k_txn[n]) #globally changes transcription rates
                
                # run simulaton (input is simulaton time)
                model.simulate(simTime,smethod='LSODA')
                
                # pull out the species from the model solution to plot
                S2 = model.output_concentration('S{2}')
                if (1+2*n > 5):
                    plt.subplot(3,5,1+2*n+4)
                    
                else:                   
                    plt.subplot(3,5,1+2*n)
                 
                plt.plot(model.t/60,(data/REP_con)*100,color='aqua',linewidth=2,linestyle='-')
                plt.plot(model.t/60,(S2/REP_con)*100,color='orange',linewidth=2,linestyle='--')
                plt.xticks(fontsize=fs)
                plt.yticks(fontsize=fs)
                plt.ylim(-10,110)
                plt.xlim(0,120)
                ax1 = plt.gca()
                ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
                #ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
                plt.xlabel('Time (min)',fontsize=fs)
                plt.ylabel('Reacted reporter (%)',fontsize=fs)
                plt.legend(['Raw Data','Calibrated'],frameon=False)
                plt.title('k='+str(k_txn[n]))
                
            
        
        







