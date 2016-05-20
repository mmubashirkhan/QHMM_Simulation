# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 12:33:37 2016

@author: mmkhan
"""

from __future__ import division

import QHMM_V6 as qm
import numpy as np
import glob as gl
import pickle
import sys
import qutip as qt
import scipy.linalg as spln
import numpy as np
import math as mt
import random as rn
import matplotlib.pyplot as plt
import sys

"""
===============================GLOBAL PARAMETERS===============================
"""
Adam = 0
Lewis = 1

if Adam == 1:
    sim_data = "Adam_sim_data"
    
if Lewis == 1:
    sim_data = "Lewis_sim_data5"

filenames_list = ['file_a','file_b','file_c','file_d',
             'file_e','file_f','file_g','file_h','file_i']

no_of_sets = 1
total_phis = 180
total_timesteps = 25

filechains_per_phi = 1000 # no of chains in a file
phis_per_file = 20

"""
===========================Count Number of Files===============================
"""
no_of_files = []
files_with_numbers = []
for f in filenames_list:
    each_file_count = 0
#    for filename in gl.glob("temp_data/**/"+str(f)+"*.pkl"):
    for filename in gl.glob(str(sim_data)+"/**/"+str(f)+"*.pkl"):
        each_file_count += 1
    files_with_numbers.append(f+" = "+str(each_file_count))
    no_of_files.append(each_file_count)
    
for f in files_with_numbers:
    print f+',000'
min_files = np.min(no_of_files)
print 'Min number of files in a set is: ', min_files

#min_files = 99
    
"""
===============================EXTRACT FILE DATA===============================
chains = ['ch0','ch1','ch2', ... ,'chn']
c = [ [ALL CHAINS PHI_1], [ALL CHAINS PHI_2] ... [ALL CHAINS PHI_n] ]
"""

chains = [['x']]*total_phis # for initialising list
file_type_count = 0

for f in filenames_list:

    print '\r     Fetching data from files ... ',f,'                \r',
    files_count = 0
    begin_phi = file_type_count*phis_per_file
    file_type_count += 1
    for filename in gl.glob(str(sim_data)+'/**/'+str(f)+'*.pkl'):
        with open(filename, 'rb') as pickle_load:
            print filename
            temp_list = pickle.load(pickle_load)
        chains = qm.extract_chains_per_file_phi(chains,temp_list,
                                        filechains_per_phi,
                                        phis_per_file,begin_phi)
        files_count += 1
        if files_count >= min_files:
            break

"""
==============================Filter Garbage [['x']]===========================
""" 
for i in range(0,total_phis):
    del chains[i][0]
"""
==============================Saving Data in Grand File========================
"""
#print 'Saving data in grand file ...'
#qm.save_list_to_file(file_name='grand_data',list_name=c)

"""
===============================Fetching Data From Grand File===================
"""
#c = qm.restore_list_from_file(dir_path='',filename_list=["grand_data.pkl"])
#c = c[0] # important
    
"""
==============First and Second Order Correlation and Gradient==================
 SOC = [ [t0...t99] [t0...t99] [t0...t99] ... [t0...t99] ] for all phis
 
 GRA = [ [P0...P179] [P0...P179] [P0...P179] ... [P0...P179] ] for all times
===============================================================================
""" 
print '\nCalculating the First and Second Order Correlation and Gradient ...'                                    

if Lewis == 1: 
    soc,grad_soc,foc,grad_foc  = qm.calc_soc_grad(chain_all = chains, total_sets=no_of_sets,
                             total_chains_per_phi = min_files*filechains_per_phi,
                             time_steps=total_timesteps, no_of_phis=total_phis)

#                     FOR ADAM'S DATA
  
if Adam == 1: 
    socAdam, gradAdam = qm.calc_soc_grad_Adam(chain_all = chains, total_sets=no_of_sets,
                                 total_chains_per_phi = min_files*filechains_per_phi,
                                 time_steps=total_timesteps, no_of_phis=total_phis)
    soc = socAdam
    grad_soc = gradAdam                             
"""
==========================Standard Deviation of FOC============================
 STDFOC = [ [t0...t99] [t0...t99] [t0...t99] ... [t0...t99] ] for all phis
===============================================================================
"""  
std_foc = []
for phi in range(0,total_phis): # for each phi
    #print phi
    f = np.array([item[phi] for item in foc])
    std_foc.append(np.std(foc[0][phi],axis=0,dtype=np.float128))    
f = []
"""
==========================Standard Deviation of SOC============================
 STDSOC = [ [t0...t99] [t0...t99] [t0...t99] ... [t0...t99] ] for all phis
===============================================================================
"""  
std_soc = []
for phi in range(0,total_phis): # for each phi
    #print phi
    s = np.array([item[phi] for item in soc])
    std_soc.append(np.std(s,axis=0,dtype=np.float128))    
s = []
"""
========================Maximum of Gradients out of many sets==================
Should be average
"""
max_gr = [] # max gradient

for time_step in range(0,total_timesteps):
    g = np.array([item[time_step] for item in grad_soc])
    max_gr.append(np.amax(np.abs(g),axis=0))
g = []

"""
========================Average of Gradients_FOC out of many sets==================
AVG_GRA = [ [P0...P179] [P0...P179] [P0...P179] ... [P0...P179] ] for all times
"""
avg_gr_foc = [] # max gradient

for time_step in range(0,total_timesteps):
    g = np.array([item[time_step] for item in grad_foc])
    avg_gr_foc.append(np.average(np.abs(g),axis=0))
g = []

"""
========================Average of Gradient_SOC out of many sets==================
AVG_GRA = [ [P0...P179] [P0...P179] [P0...P179] ... [P0...P179] ] for all times
"""
avg_gr_soc = [] # max gradient

for time_step in range(0,total_timesteps):
    g = np.array([item[time_step] for item in grad_soc])
    avg_gr_soc.append(np.average(np.abs(g),axis=0))
g = []
"""
========================Maximum of Gradients===================================
"""
max_gr_t = np.nanmax(max_gr,axis=1) # max gradient at particular timestep
max_gi = np.nanargmax(max_gr,axis=1) # corresponding index i.e. phi

"""
##==================Error Prop (FOC) for all times and max phi=======================
"""
err_prop_foc = []
for phi in range(0,total_phis):
    ep_at_phi = []
    for time_step in range(0,total_timesteps): 
        ep_at_phi.append(np.float128(std_foc[phi]) / 
                                    np.abs(np.float128(grad_foc[0][phi][time_step]))) 
    err_prop_foc.append(ep_at_phi)  
ep_at_phi = []

"""
##==================Error Prop (SOC) for all times and max phi=======================
"""
err_prop_soc = []
for phi in range(0,total_phis):
    ep_at_phi = []
    for time_step in range(0,total_timesteps): 
        ep_at_phi.append(np.float128(std_soc[phi][time_step]) / 
                                    np.float128(avg_gr_soc[time_step][phi])) 
    err_prop_soc.append(ep_at_phi)  
ep_at_phi = []

# =========================Plotting Params===================================== 
sys.exit('Testing!!!')

"""
===========================Plotting SOC and Gradient===========================
"""
for i in range(0,no_of_sets):
    print 'Data Set: ',i
    qm.plt_soc_wrt_phi(data_set=soc,set_no=i,at_timestep=2)
    qm.plt_grad_wrt_time(grad_soc,set_no=i,time_step=2)

"""
===========================Plotting FOC and Gradient===========================
"""
for i in range(0,no_of_sets):
    print 'Data Set: ',i
    qm.plt_foc_wrt_phi(data_set=foc,set_no=i,from_time=0,to_time=24)

"""
===========Plotting Error Propagation All time steps and Phi===================
"""
for i in range(0,total_phis):
    plt.plot(err_prop_foc[i])
plt.show()    
"""
===========================Plotting Error Propagation==========================
"""
for i in range(0,total_phis):
    plt.plot(err_prop_foc[i])
plt.show()

import matplotlib.pyplot as plt

phi = 120
lst2 = [item[phi] for item in avg_gr]
plt.plot(lst2,'bs',label="Gradient Avg")

plt.plot(std_soc[phi],'g^',label="Std Dev SoC")

plt.axis([1, 100, -1, 2])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()
