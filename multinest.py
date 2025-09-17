#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:09:35 2024

@author: everson
"""

import pymultinest
import json
import numpy as np
from numpy import log, exp, pi
import scipy.stats, scipy
import inspect
import os
# import mpi4py

import warnings

import tov_NoMultV2 as tov

warnings.filterwarnings('ignore', category=RuntimeWarning)

label  = 'FSU2Rsrc'
outdir = 'FSU2Rsrc'

# nuclear saturation density
rho0 = 0.15 # fm^(-3)

# likelihood - Mass
a  = 0.06 # 0.01
mu = 21


# # Set the CPU affinity for the current process
# desired_cpu_cores = [0, 1, 2]  # Set to the desired CPU cores
# os.sched_setaffinity(0, desired_cpu_cores)

def log_likelihood(cube,ndim,cnd):       
    nmemax,dmemax,rnmmax,rdmmax = tov.calculation(cube[0],1950,3.26) #frm,MDM,cnde
    loglikelihood = (nmemax+dmemax)/a - 10*mu - log(np.exp((nmemax+dmemax)/a - 10*mu) + 1)
#  			loglikelihood = (-0.5 * ((rho_trans -x0) / sigma)**(2*p))\
# 				+(-0.5 * ((R2 -12.55) / 0.57)**2) # R of M = 2.08: 12.55 +- 1.14
    return loglikelihood
 	

################# Priors ################

def prior(cube, ndim, nparams):
    # Define uniform priors within specific ranges for each parameter
    # Replace this with your own prior definition
    cube[0] = cube[0] * (0.1 - 0.005) + 0.005 # frm 
    cube[1] = cube[1] * (1960 - 1940) + 1940  # MDM
    cube[2] = cube[2] * (10.0 - 0.01) + 0.01  # cnd 

# Define the number of parameters and the parameter names
nparams=3
# ndims=1
parameter_names = ['frm','MDM','cnd']


# Set up the PyMultiNest run
pymultinest.run(log_likelihood, prior, nparams, outputfiles_basename="output/", 
                resume=False, verbose=True, init_MPI=False,n_live_points=5)


# ultranest
# ptmcmc


