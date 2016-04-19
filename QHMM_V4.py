# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:25:24 2016

@author: mmkhan
"""
from __future__ import division

import qutip as qt
import scipy.linalg as spln
import numpy as np
import math as mt
import random as rn
import matplotlib.pyplot as plt
import multiprocessing
import sys
import pickle


class QHMMv2(multiprocessing.Process):  
    
    def __init__(self):
        multiprocessing.Process.__init__(self)
        self.signal_value = ""
        self.state_vector = (qt.basis(2, 0) + qt.basis(2, 1)).unit() #unit method to again normalize the state. Operating with the number function again:
        

                
    def runGeneralQHMM(self, timeSteps,phi):

        self.current_state = self.state_vector
        self.signal_value = ''
        
        #phi = 45        
        
        self.k0 = qt.Qobj([[np.sin(phi),0],[0,np.cos(phi)]])
        self.k1 = qt.Qobj([[0,np.cos(phi)],[np.sin(phi),0]])

#        print 'k0=\n',self.k0.full
#        print 'k1=\n',self.k1.full
#        print 'K0_K0_d k1_k1_d: \n',self.k0 * self.k0.dag() + self.k1 * self.k1.dag()
#        print 'state: \n', self.state_vector
        
        for i in range(0, timeSteps):
            
            randomNo = rn.uniform(0, 1)    
            k0_psi = self.k0 * self.current_state
            k1_psi = self.k1 * self.current_state
    
            prob_0 = k0_psi.dag() * k0_psi            
            prob_0 = float(np.abs(prob_0[0])[0][0])            
            
#            print 'k0_psi \n', k0_psi.full()
#            print 'k1_psi \n', k1_psi.full()
#            print randomNo
#            print prob_0
            
            if randomNo >= np.around(prob_0, decimals=8):
                self.current_state = k0_psi.unit()
                self.signal_value += '0'
            else:
                self.current_state = k1_psi.unit()
                self.signal_value += '1'
#            print self.current_state
        return self.signal_value

#===============================Calculating probability========================
def prob_at_t(chains,t,emission_symbol):
    count = 0

    for i in range(0,len(chains)):        
        a = chains[i][t:t+len(emission_symbol)]
        b = emission_symbol

        if a == b:
            count += 1
    prob = float(count)/float(len(chains))
#    print 'count:', count, 'len: ', len(chains), 'prob:', prob
    return prob
    
def sec_ord_corr(chains,t1,emission_symbol_1,t2,emission_symbol_2):
    count = 0

    for i in range(0,len(chains)):        
        a_t1 = chains[i][t1]
        b_t1 = emission_symbol_1
        
        a_t2 = chains[i][t2]
        b_t2 = emission_symbol_2

        if a_t1 == b_t1 and a_t2 == b_t2:
            count += 1
    return (float(count)/float(len(chains)))
#============================saving data to a file=============================
def save_list_to_file(file_name,list_name):
    with open(file_name+'.pkl', 'wb') as pickle_file:
        pickle.dump(list_name, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
    return 0
#============================Extract Data======================================        
def restore_list_from_file(file_name):
    with open(file_name+'.pkl', 'rb') as pickle_load:
        new_list = pickle.load(pickle_load)
    return new_list  
#========================Reformatting/Gathering Data===========================
def extract_chains_per_phi(chain_all):
    new_chain_all = []
    for chains in chain_all:
        new_chain = []    
        for c in chains:
            for b in c:
                new_chain.append(b)
        new_chain_all.append(new_chain)
    return new_chain_all
#==============================================================================
phi_0 = []
phi_1 = []
phi_all = []
grad_2_corr = []
chain_all = []
#==========================Set Coefficients=For all Chains=====================

file_name = 'file_b'
#angles = np.linspace(0.5,39.5,40)  # file_a
angles = np.linspace(40.5,79.5,40)  # file_b
#angles = np.linspace(80.5,119.5,40)  # file_c
#angles = np.linspace(120.5,159.5,40)  # file_d
#angles = np.linspace(160.5,199.5,40)   # file_e
#angles = np.linspace(200.5,339.5,40)  # file_f
#angles = np.linspace(240.5,279.5,40)  # file_g
#angles = np.linspace(280.5,319.5,40)  # file_h
#angles = np.linspace(320.5,359.5,40)  # file_i

x = []

for i in range(0,500):
    x.append(QHMMv2())
    
for i in angles:
    print i
    phi_all.append(i)

    time_steps = 100
    no_of_chains = 10
    
    chains = []
    for j in range(0,no_of_chains):
        print j

        for sim in x:
            chains.append(sim.runGeneralQHMM(time_steps,np.deg2rad(i)))
        
    for c in chains:
        chain_all.append(c)

#==============================================================================
save_list_to_file(file_name,chain_all)

lst_a = restore_list_from_file(file_name)
#new_chain_per_phi = []
#new_chain_per_phi = extract_chains_per_phi(lst_a)    

#sys.exit('Stop!!!')
##==================================1st sim===================================== 

##    
#first_order_corr_all = []
#second_order_corr_all = []
# 
#for i in range(0,np.size(angles)):
#    first_order_corr = []
#    second_order_corr = [] 
#    for k in range(0,time_steps):
#        foc = prob_at_t(chain_all[i],k,'1')
#        first_order_corr.append(foc)
#        second_order_corr.append(sec_ord_corr(chain_all[i],k,'1',time_steps-1,'1'))
#
#    first_order_corr_all.append(first_order_corr)
#    second_order_corr_all.append(second_order_corr) 
#    
###========================STANDARD DEVIATION====================================
##for every time step we calculate Std Dev
#    
#std_dev_all = []
#for i in range(0,time_steps):
#    lst2 = [item[i] for item in second_order_corr_all]
#    std_dev_all.append(np.std(lst2))
#
###======================== Gradient of SOC wrt Phi for each time step===========
#gradient_all = []
#diff_all = []
#for i in range(0,time_steps):
#    lst2 = [item[i] for item in second_order_corr_all]
#    gradient_all.append(np.gradient(lst2,np.size(angles)))
#    diff_all.append(np.diff(lst2))
#===============================================================================    
#del_phi_all = []
#for i in range(0,no_of_chains):
#    lst2 = [item[i] for item in second_order_corr_all]
#    del_phi = []
#    for j in range(0,np.size(lst2)):
#        del_phi.append(std_deviation[i]/np.abs(lst2[j]))
#    del_phi_all.append(del_phi)
#
#col_del_phi = []
#for i in range(0,np.size(del_phi_all[0])):
#    col_del_phi.append([item[i] for item in del_phi_all])
#
#plt.plot(col_del_phi)
