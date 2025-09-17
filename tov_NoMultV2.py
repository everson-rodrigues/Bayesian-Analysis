#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This code uses bisection to calculate a TOV. Some information is 
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
###
###
###
#### Part of the code omitted for confidentiality reasons.
###
###
###
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


###
###
###
#### Part of the code omitted for confidentiality reasons.
###
###
###
        
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


