#!/usr/bin/env python
'''
STU analisis for SDFDM model arXiv:1411.1335
TODO: Only T included
return: array of T paremeter
'''
from numpy import *
import commands
import numpy as np
import sys
import scipy as sp
import scipy.integrate as integrate
from scipy.integrate  import quad
import scipy.optimize

## Calculate B0 for p**2=0
def B00(m1,m2):
    bb00=1.0-(m1**2*log(m1**2)-m2**2*log(m2**2))/(m1**2-m2**2)
    return bb00

# Function to calculate B0
def logDelta(x,pp,m1,m2):
    yy=x*m1**2+(1.0-x)*m2**2-x*(1.0-x)*pp**2
    return np.log(yy)

## Calculate B0 for p**2 neq 0
def B0(pp,m1,m2):
    yy=-np.asarray(quad( logDelta, 0.0, 1.0, args=(pp,m1,m2) ))
    return yy[0]

# Function to calculate H
def f1(x,pp,m1,m2):
    f2=( 4.0*x*(1.0-x)*pp**2 - 2.0*(1.0-x)*m1**2 - 2.0*x*m2**2)*logDelta(x,pp,m1,m2)
    return f2

## Calculate H for arbitrary p**2
def HH(pp,m1,m2):
    hh=-np.asarray(quad( f1, 0.0, 1.0, args=(pp,m1,m2) ))
    return hh[0]



def stu(ms,md,ll,tbeta,vev=246,sw=sqrt(0.23),mz=91.18,gg=0.653):
            
    r2=np.sqrt(2.0)
    cw=sqrt(1.-sw**2)
    mw=mz*cw
    ee=gg*sw
    alpha=ee**2/(4.0*np.pi)

    yy=ll*cos(arctan(tbeta))  #random.uniform(0.1,10)
    yyp=ll*sin(arctan(tbeta))# random.uniform(0.1,10)
 
    #Mass matrix from arXiv:1411.1335, eq(12). Fix a global sign in Lagrangian, e.g L -> -L and MN -> -ms
    Mmass=[[-ms, -yy*vev/r2, yyp*vev/r2],[yy*vev/r2, 0.0, -md], [-yyp*vev/r2, -md, 0]]
    eigsys=np.linalg.eigh(Mmass)
    eigva=eigsys[0]
    UU0=eigsys[1]
 
    if(abs(eigva[0])<abs(eigva[1])and abs(eigva[0])<abs(eigva[2])):
        masses=eigva
        UU=UU0
    if(abs(eigva[1])<abs(eigva[0])and abs(eigva[1])<abs(eigva[2])):
        vecs=np.transpose(UU0)
        VV=[vecs[1], vecs[0], vecs[2]]
        UU=np.transpose(VV)
        masses=[eigva[1],eigva[0],eigva[2]]
    if(abs(eigva[2])<abs(eigva[0])and abs(eigva[2])<abs(eigva[1])):
        vecs=np.transpose(UU0)
        VV=[vecs[2], vecs[0], vecs[1]]
        UU=np.transpose(VV)
        masses=[eigva[2],eigva[0],eigva[1]]
 
    #Inputs for the ST parameters: eq(16) from arXiv:1411.1335. 
    Cli = np.asarray([UU[2][0], UU[2][1], UU[2][2]])/r2
    Cri = -1.0*np.asarray([UU[1][0], UU[1][1], UU[1][2]])/r2
    Nlij=[[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
    for ii in range(3):
        for jj in range(3):
            Nlij[ii][jj]=-0.5*(UU[2][ii]*UU[2][jj]-UU[1][ii]*UU[1][jj])
    Nrij=-1.0*np.asarray(Nlij)
 
    #Calculating T parameter
    #W vacuum polarization 
    PiWW0=0
    chk=0
    for ii in range(3):
        PiWW0=PiWW0-gg**2/(16.0*np.pi**2)*( (Cli[ii]**2+Cri[ii]**2)*HH(0.0,masses[ii],md) + 4.0*Cli[ii]*Cri[ii]*md*masses[ii]*B0(0.0,masses[ii],md) )
        chk=chk+Cli[ii]*Cri[ii]*md
    
    #Z vacuum polarization 
    PiZZ01=-gg**2*(1.0-2.0*sw**2)**2/(32.0*np.pi**2*cw**2)*( HH(0.0,md,md) + 2.0*md**2*B0(0.0,md,md))
 
    PiZZ02=0
    for ii in range(3):
            for jj in range(3):
                PiZZ02=PiZZ02-gg**2/(32.0*np.pi**2*cw**2)*( ( Nlij[ii][jj]**2 + Nrij[ii][jj]**2 )*HH(0.0,masses[ii],masses[jj]) )
                
    PiZZ03=0
    for ii in range(3):
        for jj in range(3):
            PiZZ03=PiZZ03-gg**2/(32.0*np.pi**2*cw**2)*4.0*( Nlij[ii][jj]*Nrij[jj][ii]*masses[ii]*masses[jj]*B0(0.0,masses[ii],masses[jj]) )
 
    #T parameter        
    TT=1.0/(alpha*mw**2)*( PiWW0 - (PiZZ01+PiZZ02+PiZZ03)*cw**2 )
   
    return TT




stumap=np.vectorize(stu,excluded={'vev':246.2,'sw':sqrt(0.23),'mz':91.18,'gg':0.653}\
                    ,doc='Input for pyfunc below: MN,MDF,lambda,tanb,vev=266.2,sw=sqrt(0.23),mz=91.18,gg=0.653')

#####-------------------------####
'''
STU analisis for SDFDM model  arXiv:0705.4493
TODO: Only T included
return: array of T paremeter
'''

def PI0(m1,m2,LAMBDA=1e+16):
    if m1==m2:
        pizero=0.
    elif m1==-m2:
        pizero=1.0/(16.0*np.pi**2)*( 4.0*m1**2*(-1+np.log(LAMBDA**4/m1**4)) )
    else:
        pizero = 1./(16*np.pi**2)*( (m1-m2)**2*np.log(LAMBDA**4/(m1**2*m2**2))-2*m1*m2\
            +( 2*m1*m2*(m1**2+m2**2)-m1**4-m2**4 )/( m1**2-m2**2 )*np.log(m1**2/m2**2) )
    return pizero

def A(m1,m2,v=246.2,alpha_em=1./128.,LAMBDA=1e+16,):
    aa=1./(alpha_em*v**2)*PI0(m1,m2,LAMBDA)   
    return 2*aa  #Note the factor 2.0 because we are dealing with Majorana fermions

def stu_new(ms,md,ll,tbeta,vev=246.2,sw=np.sqrt(0.23),gg=0.653,LAMBDA=1E16):
    alpha=(gg*sw)**2/(4.0*np.pi)
    yy=ll*np.cos(np.arctan(tbeta))  
    yyp=ll*np.sin(np.arctan(tbeta)) 
    gamma=(yyp+yy)/2.0
    beta=(yyp-yy)/2.0
    
    #Mass matrix from arXiv:1411.1335, eq(12). 
    Mmass=[[-md,              0,               -beta*vev ],
           [0.0,              md,              -gamma*vev], 
           [-beta*vev, -gamma*vev, ms]]
    #Mass matrix from PRD by ERAMO, eq(6) but with md->-md, ms->-2ms in 
    #order to be compatible with arXiv:1411.1335, eq(12).   
    #These changes were also implemented in the T parameter formulae.

    eigsys=np.linalg.eigh(Mmass)
    m=eigsys[0]
    V=eigsys[1]
    
    TT1=0
    TT2=0
    for i in range(3):
        TT1=TT1+(V[0,i])**2*A(-md,m[i],vev,alpha,LAMBDA)+(V[1,i])**2*A(-md,-m[i],vev,alpha,LAMBDA)
        for j in range(3):
            TT2=TT2+(V[0,i]*V[1,j]+V[1,i]*V[0,j])**2*A(m[i],-m[j],vev,alpha,LAMBDA)

    TT=TT1-0.5*TT2

    return TT

stumap_new=np.vectorize(stu_new,excluded={'vev':246.2,'sw':np.sqrt(0.23),'gg':0.653,'LAMBDA':1E16},\
                      doc='Input for pyfunc below: MN,MDF,lambda,tanb,vev=266.2,sw=sqrt(0.23),gg=0.653')
