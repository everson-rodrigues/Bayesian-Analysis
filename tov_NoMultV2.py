#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 14:58:31 2024

@author: everson
"""

import numpy as np
import pandas as pnd
import scipy.integrate as spi
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import csv
import gc
import eos

global data, eNM, pNM, eDM, pDM, eNM_datastars, eDM_datastars, pNM_datastars, pDM_datastars
global pn, pd, mn, md, nme, dme, rnm, rdm, nme, dme, frc, frm
global search, search_fr

start_time = time.time()
# Constants
gms = 1.475562  # [L]/[M]
ms = 5660.57
h = 197.3269631

nme_frm = []
dme_frm = []
rnm_frm = []
rdm_frm = []
totalmass_frm = []
totalradius_frm = []
starradius_frm = []
list_fr = []
index_fr = []
deltaradius_frm = []

nstars = 8000
nsteps = 1
niterations = 20
search_fr = 1e-4

r_start = 1e-6
r_end = 100
step = 1
search = 5e-9

rlist = [[] for i in range(nstars)]
rnm = np.zeros([nstars])
rdm = np.zeros([nstars])
nme = np.zeros([nstars])
dme = np.zeros([nstars])
frc = np.zeros([nstars])

num_processes = 1  # Changed to 1 for no parallelization

manager = None

pn = [0] * nstars
pd = [0] * nstars
mn = [0] * nstars
md = [0] * nstars
rlist = [0] * nstars

nme = [0] * nstars
dme = [0] * nstars
rnm = [0] * nstars
rdm = [0] * nstars
frc = [0] * nstars


def compute_tov(i, pn, pd, mn, md, nme, dme, rnm, rdm, rlist,
                pNM_datastars, pDM_datastars, efpNM, efpDM):
    init = np.array([pNM_datastars, pDM_datastars, 0, 0])
    sol = spi.solve_ivp(lambda r, y: tov_two_fluids(r, y, efpNM, efpDM), [r_start, r_end], init,
                        method='RK23', max_step=step, events=event,
                        atol=1e-6, rtol=1e-6)

    rlist[i] = sol.t
    pn[i] = sol.y[0]
    pd[i] = sol.y[1]
    mn[i] = sol.y[2]
    md[i] = sol.y[3]

    if min(pn[i]) < search:
        rnm[i] = rlist[i][np.where(pn[i] < search)[0][0]]
        nme[i] = mn[i][np.where(pn[i] < search)[0][0]]

    if min(pd[i]) < search:
        rdm[i] = rlist[i][np.where(pd[i] < search)[0][0]]
        dme[i] = md[i][np.where(pd[i] < search)[0][0]]

    else:
        rnm[i] = np.where(pn[i] > min(pn[i]))[0][-1]
        rdm[i] = np.where(pd[i] > min(pd[i]))[0][-1]


# TOV equations for two fluids
def tov_two_fluids(r, y, efpNM, efpDM):
    # m = total mass, pn = pressure for normal matter, pd = pressure for dark matter
    pn, pd, mn, md = y
    
    dpndr = - (pn+efpNM(pn)) * ((mn+md) + 4 * np.pi *
                                r**3 * (pn+pd)/ms)/(r**2/gms-2*r*(mn+md))
    dpddr = - (pd+efpDM(pd)) * ((mn+md) + 4 * np.pi *
                                r**3 * (pn+pd)/ms)/(r**2/gms-2*r*(mn+md))
    dmndr = 4 * np.pi * r**2 * (efpNM(pn)/ms)
    dmddr = 4 * np.pi * r**2 * (efpDM(pd)/ms)

    return [dpndr, dpddr, dmndr, dmddr]


# Event to stop the solving of equations
def event(r, y):
    if abs(y[0]) > abs(y[1]):
        ev = y[0]
    else:
        ev = y[1]
    return ev-search/2


def findminimum(pn, pd, mn, md, nme, dme, rnm, rdm, rlist, 
                pNM_datastars, pDM_datastars, efpNM,efpDM):
    results = []
    ni = int(0.1*nstars) #100 #5
    nf = nstars-1
    nm = int((nf+ni)/2.)
    # ni = nm = nf = 150
    nim = int((nm+ni)/2.)
    nmf = int((nm+nf)/2.)
    nit = 20
    nsteps = 1
    while nsteps <= nit:
        # print(f" Finding n_min: {nsteps}", end="")
        for i in [ni, nim, nm, nmf, nf]:
            # if frc[i] == 0:
            compute_tov(i, pn, pd, mn, md, nme, dme, rnm, rdm, rlist, 
                        pNM_datastars[i], pDM_datastars[i], efpNM, efpDM)
 
        # Access the results
        for i in [ni, nim, nm, nmf, nf]:
            frc[i] = dme[i]/(nme[i]+dme[i])
        # print("It, FR, Radius, Mass")
        # print(ni, '{:.5f}'.format(frc[ni]), '{:.5f}'.format(rnm[ni]) if rnm[ni] > rdm[ni] else '{:.5f}'.format(rnm[ni]),'{:.5f}'.format( nme[ni]+dme[ni]))
        # print(nim, '{:.5f}'.format(frc[nim]), '{:.5f}'.format(rnm[nim]) if rnm[nim] > rdm[nim] else '{:.5f}'.format(rnm[nim]),'{:.5f}'.format( nme[nim]+dme[nim]))
        # print(nm, '{:.5f}'.format(frc[nm]), '{:.5f}'.format(rnm[nm]) if rnm[nm] > rdm[nm] else '{:.5f}'.format(rnm[nm]), '{:.5f}'.format(nme[nm]+dme[nm]))
        # print(nmf, '{:.5f}'.format(frc[nmf]), '{:.5f}'.format(rnm[nmf]) if rnm[nmf] > rdm[nmf] else '{:.5f}'.format(rnm[nmf]),'{:.5f}'.format( nme[nmf]+dme[nmf]))
        # print(nf, '{:.5f}'.format(frc[nf]), '{:.5f}'.format(rnm[nf]) if rnm[nf] > rdm[nf] else '{:.5f}'.format(rnm[nf]), '{:.5f}'.format(nme[nf]+dme[nf]))
        # print(abs(frm-frc[ni]), abs(frm-frc[nm]), abs(frm-frc[nf]))
        
        if frc[nim] < frc[ni] and frc[nim] < frc[nm]:
            ni_new = ni
            nf_new = nm
            nm_new = int((ni_new+nf_new)/2)
            nim = int((ni_new+nm_new)/2.)
            nmf = int((nm_new+nf_new)/2.)
        elif frc[nmf] < frc[nm] and frc[nmf] < frc[nf]:
            ni_new = nm
            nf_new = nf
            nm_new = int((ni_new+nf_new)/2)
            nim = int((ni_new+nm_new)/2.)
            nmf = int((nm_new+nf_new)/2.)
        else:
            if frc[nim]<frc[nmf]:
                ni_new = ni
                nf_new = nm
                nm_new = int((ni_new+nf_new)/2)
                nim = int((ni_new+nm_new)/2.)
                nmf = int((nm_new+nf_new)/2.)
            elif frc[nmf]<frc[nim]:
                ni_new = nm
                nf_new = nf
                nm_new = int((ni_new+nf_new)/2)
                nim = int((ni_new+nm_new)/2.)
                nmf = int((nm_new+nf_new)/2.)
                
        ni = ni_new
        nm = nm_new
        nf = nf_new
        nsteps = nsteps+1
        if abs(frc[nf]-frc[nm])<=search_fr: break
    
    if frc[nm]<frc[nf] and frc[nm]<frc[ni]: nmin=nm
    if frc[ni]<frc[nf] and frc[ni]<frc[nm]: nmin=ni
    if frc[nf]<frc[nm] and frc[nf]<frc[ni]: nmin=nf
    print(" ")  
    print ('Minimum - frc=','{:.5f}'.format(frc[nm]),'  nmin=', nmin)
    return frc[nm],nmin


def calculation(frm,MDM,cnde):
    # print("--------------------------- start ------------------------------------------")
    # print (frm,MDM,cnde)
    n0=0.15
    fr=1.4
    frf=1.5
    csde=0
    SRe=1
    # Import the data for the equation of state    
    while fr<frf:
            # print(" ")  
            # print(" ")  
            # print("------------------ FR --------------------")
            # print(fr)
            nme = [0] * nstars
            dme = [0] * nstars
            rnm = [0] * nstars
            rdm = [0] * nstars
            frc = [0] * nstars
        
            frstr=str("{:.3f}".format(fr))
            # print("--------- EOS ---------")
            data = eos.EoScalculation(fr, MDM, cnde, csde, SRe, 0.001*n0)
            data = np.transpose(np.transpose(data)[5:-3])
            
            data=pnd.DataFrame(np.transpose(data))
            data=data[::-1]
            rho = data.iloc[:, 0]
            
            eNM = data.iloc[:, 1]
            pNM = data.iloc[:, 2]
            
            eDM = data.iloc[:, 3]
            pDM = abs(data.iloc[:, 4])
            
            
            # Interpolate the pressure as function of the energy densitye and vice-versa
            efpNM=interpolate.interp1d(pNM,eNM,fill_value='extrapolate', kind='linear')
            efpDM=interpolate.interp1d(pDM,eDM,fill_value='extrapolate', kind='linear')
            
            pfeNM=interpolate.interp1d(eNM,pNM,fill_value='extrapolate', kind='linear')
            pfeDM=interpolate.interp1d(eDM,pDM,fill_value='extrapolate', kind='linear')
            
            # Data Stars 
            # eNM_datastars=np.exp(np.linspace(np.log(eNM[1]), eNM[eNM.shape[0]-1], nstars))
            # eDM_datastars=np.exp(np.linspace(np.log(eDM[1]), eDM[eDM.shape[0]-1], nstars))
            eNM_datastars=np.linspace(eNM[1], eNM[eNM.shape[0]-1], nstars)
            eDM_datastars=np.linspace(eDM[1], eDM[eDM.shape[0]-1], nstars)
            
            pNM_datastars=pfeNM(eNM_datastars)
            pDM_datastars=pfeDM(eDM_datastars)
        
    
            nmin=51
         
            
            def store(nme, dme, rnm, rdm, nres):
                # starradius = np.where(rnm >= rdm, rnm, rdm)
                # totalmass = nme+dme
                nme_frm.append(nme[nres])
                dme_frm.append(dme[nres])
                rnm_frm.append(rnm[nres])
                rdm_frm.append(rdm[nres])
                totalmass_frm.append(nme[nres]+dme[nres])
                starradius_frm.append(rdm[nres] if rdm[nres] > rnm[nres] else rnm[nres])
                deltaradius_frm.append(abs(rnm[nres]-rdm[nres]))
                
                return totalmass_frm, starradius_frm,nme_frm, rnm_frm, dme_frm, rdm_frm, deltaradius_frm
            
            # Necessary to stop the integration given the event
            event.terminal = True
            
            # range of radius of the star
            # print("--------- Finding Minimun ---------")    
            [frcmin,nmin]=findminimum(pn, pd, mn, md, nme, dme, rnm, rdm, rlist, 
                                      pNM_datastars, pDM_datastars, efpNM,efpDM)
            nsteps = 1
            cont=0
            
            ni = int(0.1*nstars)#100 #5
            nf = nmin
            nm = int(abs(nf-ni)/2)
            
            while nsteps <= niterations:
                for i in [ni, nm, nf]:
                    if frc[i] == 0:
                        compute_tov(i, pn, pd, mn, md, nme, dme, rnm, rdm, rlist,
                                    pNM_datastars[i], pDM_datastars[i], efpNM,efpDM)
 
                # Access the results
                for i in [ni, nm, nf]:
                    frc[i] = dme[i]/(nme[i]+dme[i])
        
                if frm > frc[int(0.1*nstars)] and frm > frc[nstars-1] and nsteps == 1:
                    break_var=1
                    break
                if frm < frcmin:
                    break_var=1
                    break
             
                list_fr.append(round(fr, 2))
                # index_fr.append(nres) 
                
                if cont==0:     
                    if frm < frc[nm]:
                        nmnew = int((nf+nm)/2)
                        ninew = nm
                        nfnew = nf
                    if frm > frc[nm]:
                        nmnew = int((ni+nm)/2)
                        ninew = ni
                        nfnew = nm            
                if cont==2:     
                    if frm < frc[nm]:
                        nmnew = int((ni+nm)/2)
                        ninew = ni
                        nfnew = nm
                    if frm > frc[nm]:
                        nmnew = int((nm+nf)/2)
                        ninew = nm
                        nfnew = nf
                        
                # print('------------')         
                # print("It, FR, Radius, Mass")
                # print(ni, frc[ni], rnm[ni] if rnm[ni] > rnm[ni] else rnm[ni], nme[ni]+dme[ni])
                # print(nm, frc[nm], rnm[nm] if rnm[nm] > rnm[nm] else rnm[nm], nme[nm]+dme[nm])
                # print(nf, frc[nf], rnm[nf] if rnm[nf] > rnm[nf] else rnm[nf], nme[nf]+dme[nf])
                # print(abs(frm-frc[ni]), abs(frm-frc[nm]), abs(frm-frc[nf]))
        
                if abs(frm-frc[ni]) < search_fr:
                    nres = ni
                    list_fr.append(round(fr, 2))
                    # index_fr.append(nres)
                    # store(nme, dme, rnm, rdm, nres)
                    store(nme, dme, rnm, rdm, ni)
                    store(nme, dme, rnm, rdm, nm)
                    store(nme, dme, rnm, rdm, nf)
                    mt=nme[nm]+dme[nm]
                    print ('--------------------------------------------' '\n' 
                           'frm=',frm,'MDM=',MDM,'cnde=',cnde,' nmin=', nmin, 'frmin=', '{:.5f}'.format(frc[nmin]), '\n' 
                           'nm=',nm,'frc=', '{:.5f}'.format(frc[nm]),'r=','{:.5f}'.format(rnm[nm]) if rnm[nm] > rnm[nm] else '{:.5f}'.format(rnm[nm]),
                           'm=',  '{:.5f}'.format(mt),'\n')
                    cont=cont+1
                    # print(nres)
    
    
                elif abs(frm-frc[nm]) < search_fr:
                    nres = nm
                    list_fr.append(round(fr, 2))
                    # index_fr.append(nres)
                    # store(nme, dme, rnm, rdm, nres)
                    # store(nme, dme, rnm, rdm, nres)
                    store(nme, dme, rnm, rdm, ni)
                    store(nme, dme, rnm, rdm, nm)
                    store(nme, dme, rnm, rdm, nf)    
                    mt=nme[nm]+dme[nm]
                    print ('--------------------------------------------' '\n' 
                           'frm=',frm,'MDM=',MDM,'cnde=',cnde,' nmin=', nmin, 'frmin=', '{:.5f}'.format(frc[nmin]), '\n' 
                           'nm=',nm,'frc=', '{:.5f}'.format(frc[nm]),'r=','{:.5f}'.format(rnm[nm]) if rnm[nm] > rnm[nm] else '{:.5f}'.format(rnm[nm]),
                           'm=',  '{:.5f}'.format(mt),'\n')
                    cont=cont+1
                    # print(nres)
    
                elif abs(frm-frc[nf]) < search_fr:
                    nres = nf
                    list_fr.append(round(fr, 2)) 
                    # index_fr.append(nres)
                    # store(nme, dme, rnm, rdm, nres)
                    store(nme, dme, rnm, rdm, ni)
                    store(nme, dme, rnm, rdm, nm)
                    store(nme, dme, rnm, rdm, nf) 
                    mt=nme[nm]+dme[nm]
                    print ('--------------------------------------------' '\n' 
                           'frm=',frm,'MDM=',MDM,'cnde=',cnde,' nmin=', nmin, 'frmin=', '{:.5f}'.format(frc[nmin]), '\n' 
                           'nm=',nm,'frc=', '{:.5f}'.format(frc[nm]),'r=','{:.5f}'.format(rnm[nm]) if rnm[nm] > rnm[nm] else '{:.5f}'.format(rnm[nm]),
                           'm=',  '{:.5f}'.format(mt),'\n')
            
                    cont=cont+1
                else: nres=1
    
                    # print(nres)
        
                ni = ninew
                nm = nmnew
                nf = nfnew
                nsteps = nsteps+1
                
                if cont==1:
                    cont=cont+1
                    ni=nmin
                    nf=nstars-1
                    nm=int((nf-ni)/2)
                    if frm>frc[nf]: break
                    print("------------------------------------------------------------------------------")    
                    
                if cont>=3:
                    break
                


                # print ('frm=',frm,'MDM=',MDM,'cnde=',cnde,'fr=',fr,'frc=','{:.5f}'.format(frc[nm]),' nmin=', nm)
            # Convert manager lists to regular lists
            # pn = list(pn)
            # pd = list(pd)
            # mn = list(mn)
            # md = list(md)
            # rlist = list(rlist)
        
            # nme = np.array(nme)
            # dme = np.array(dme)
            # rnm = np.array(rnm)
            # rdm = np.array(rdm)
        
            end_time = time.time()
            elapsed_time = end_time - start_time
            
            
            # if fr<=0.96*frf:
            #     fr=fr+0.02
            # elif fr>0.96*frf and fr<0.98*frf:
            #     fr=fr+0.04
            # elif fr>=0.98*frf:
            fr=fr+0.04
    
            exportedlist = [totalmass_frm, starradius_frm,
                        nme_frm, rnm_frm, dme_frm, rdm_frm, 
                        deltaradius_frm]
            
            if 'break_var' in locals(): 
                if break_var==1:
                    nme[0],dme[0],rnm[0],rdm[0],nres = 0,0,0,0,0
                    store(nme, dme, rnm, rdm, nres)
                
            nmemax = max(nme_frm)
            dmemax = max(dme_frm)
            rnmmax = max(dme_frm)
            rdmmax = max(rdm_frm)

    # print("--------------------------- end ------------------------------------------")
    return nmemax,dmemax,rnmmax,rdmmax
        # Your list of data ["MT", "R", "MNM", "RNM", "MDM", "RDM", "RDM-RNM"]
        
# calculation(0.02,1950,3.26)     


