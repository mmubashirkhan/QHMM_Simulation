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


class QHMM(multiprocessing.Process):  
    
    def __init__(self):
        multiprocessing.Process.__init__(self)
        self.signal_value = ""
        self.state_vector = (qt.basis(2, 0) + qt.basis(2, 1)).unit() #unit method to again normalize the state. Operating with the number function again:
                
    def runGeneralQHMM(self, timeSteps,phi):

        self.current_state = self.state_vector
        self.signal_value = ''       
        
# Kraus Operator Adam
#        self.k0 = qt.Qobj([[np.cos(phi),0],[0,np.cos(phi)]])
#        self.k1 = qt.Qobj([[0,np.sin(phi)],[np.sin(phi),0]])

# Kraus Operator Lewis
#        self.k0 = qt.Qobj([[np.cos(phi),0],[0,np.sin(phi)]])
#        self.k1 = qt.Qobj([[0,np.cos(phi)],[np.sin(phi),0]])


# Kraus Operator with Gamma
#        del_t = 10**-4
#        Note:- Here gamma is written as phi  
#        self.k0 = qt.Qobj([[1,0],[0, np.exp( (-1/2) * phi * del_t) ]])
#        self.k1 = qt.Qobj([[0,0],[0,np.sqrt( phi * del_t ) ]])

#Kraus Operator
        del_t = 10**-4
        gamma = 1
#       Note:- Here gamma is taken   
        self.k0 = qt.Qobj([[1,0],[0, np.exp( (-1/2) * gamma * del_t) ]])
        self.k1 = qt.Qobj([[0, np.sqrt( gamma * del_t )*np.sin(phi) ],
                           [0, np.sqrt( gamma * del_t )*np.cos(phi) ]])
        
        for i in range(0, timeSteps):
            
            randomNo = rn.uniform(0, 1)    
            k0_psi = self.k0 * self.current_state
            k1_psi = self.k1 * self.current_state
    
            prob_0 = k0_psi.dag() * k0_psi            
            prob_0 = float(np.abs(prob_0[0])[0][0])            
            
            if randomNo <= np.around(prob_0, decimals=8):
                self.current_state = k0_psi.unit()
                self.signal_value += '0'
            else:
                self.current_state = k1_psi.unit()
                self.signal_value += '1'
        return self.signal_value
#===============================GENERATE DATA==================================
def generateChains(phi_start,phi_end,phi_total,time_steps,no_of_instances,chains_per_instance):

    angles = np.linspace(phi_start,phi_end,phi_total)

    phi_all = []
    chain_all = []    
    x = []
    
    for i in range(0,no_of_instances):
        x.append(QHMM())
        
    for i in angles:
        print '\r','phi=',i+0.5, '           \r',

        phi_all.append(i)
    
        #time_steps = 100
        no_of_chains = chains_per_instance
        
        chains = []
        for j in range(0,no_of_chains):    
            for sim in x:
                chains.append(sim.runGeneralQHMM(time_steps,np.deg2rad(i)))
            
        for c in chains:
            chain_all.append(c)
    return chain_all

#===============================Calculating probability========================
def prob_at_t(chains,t,emission_symbol):
    count = 0

    for i in range(0,len(chains)):        
        a = chains[i][t]
        
        b = emission_symbol

        if a == b:
            count += 1
    prob = float(count)/float(len(chains))
    return prob
    
def sec_ord_corr(chains,t1,emission_symbol_1,t2,emission_symbol_2):
    
    ones_at_t1_t2 = 0
    ones_at_t1 = 0
    ones_at_t2 = 0
    
    for i in range(0,len(chains)): 
        a_t1 = chains[i][t1]
        b_t1 = emission_symbol_1
        
        a_t2 = chains[i][t2]
        b_t2 = emission_symbol_2

        if a_t1 == b_t1:
            ones_at_t1 += 1
        
        if a_t2 == b_t2:
            ones_at_t2 += 1

        if a_t1 == b_t1 and a_t2 == b_t2:
            ones_at_t1_t2 += 1
    try:
        soc = (np.float128(ones_at_t1_t2) / ones_at_t1) /  (np.float128(ones_at_t2) / len(chains))
    except ZeroDivisionError:
        soc = np.Inf
    return soc

def sec_ord_corr_Adam(chains,t1,emission_symbol_1,t2,emission_symbol_2):
    
    ones_at_t1_t2 = 0
    
    for i in range(0,len(chains)): 
        a_t1 = chains[i][t1]
        b_t1 = emission_symbol_1
        
        a_t2 = chains[i][t2]
        b_t2 = emission_symbol_2

        if a_t1 == b_t1 and a_t2 == b_t2:
            ones_at_t1_t2 += 1
    soc = np.float128(ones_at_t1_t2) / len(chains)
    
#    print '\r     ',np.float128(ones_at_t1_t2),'--',len(chains),'                                           \r',
#    sys.exit('soc test!!!')
    return soc
#============================saving data to a file=============================
def save_list_to_file(file_name,list_name):
    with open(file_name+'.pkl', 'wb') as pickle_file:
        pickle.dump(list_name, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
    return 0
#============================Extract Data======================================        
def restore_list_from_file(dir_path,filename_list):
    new_list = []
    for file_name in filename_list:
        print '\rFetching data from files ... ',file_name,'                                           \r',
        with open(dir_path+file_name, 'rb') as pickle_load:
            temp_list = pickle.load(pickle_load)
        new_list.append(temp_list)
    return new_list  
#========================Reformatting/Gathering Data===========================
def extract_chains_per_phi(lst,chains_per_angle):
    lst_a = []
    total_phis = 180
    a = 0
    d = chains_per_angle
    for n in range(1,total_phis+1):
        term = a + (n-1)*d         
        lst_a.append(lst[term:term+chains_per_angle])
    return lst_a
    
def extract_chains_per_file_phi(chain_all,file_lst,chains_per_phi,phis_per_file,begin_phi):
    for n in range(0,phis_per_file):
        term = n*chains_per_phi
        chain_all[n+begin_phi] = chain_all[n+begin_phi] + file_lst[term:term + chains_per_phi]
    return chain_all
#==============================================================================
def fetch_data_from_files(chains_per_set,dir_path,file_name_list):# per Phi chains_per_set = 250 for 10sets or 2500 for one set
    lst_all = restore_list_from_file(dir_path,file_name_list)
    lst = []
    for i in lst_all:
        for j in i:    
            lst.append(j)
    return extract_chains_per_phi(lst,chains_per_set)    

""" 
=====================Calc Second Order Corr and GRADIENT=======================
This function takes chain_all as all chains for all phi's and sets
chain_all = [ [all chains for Phi_0] [all chains for Phi_1] ... [all chains for Phi_n] ]
It generates sets of equal sizes and returns those sets with SOC and GRAD
"""
def calc_soc_grad(chain_all,total_sets,total_chains_per_phi,time_steps,no_of_phis):
    print '\r Calculating Gradient and Second Order Correlations ...                                                   \r',
    chains_per_phi_per_set = int(total_chains_per_phi / total_sets) 
    print 'chains_per_phi_per_set = ', chains_per_phi_per_set

    fixed_t1 = 1

    soc_for_all_sets = []
    foc_for_all_sets = []
    gradsoc_for_all_sets = []
    gradfoc_for_all_sets = []

    for step in range(0,total_chains_per_phi,chains_per_phi_per_set): # for every sub-set
        print '\r                   Processing Set No starting from:',step,'                                           \r',
        soc_for_all_timesteps = []
        foc_for_all_timesteps = []
        foc_for_all_phis = []
        for i in range(0,no_of_phis): # for every phi
        ###========================SECOND ORDER CORR===========================   
            print '\r',i, '          \r',      
            second_order_corr = []
            first_order_corr = []
            
            for t2 in range(0,time_steps):
                soc = sec_ord_corr(chain_all[i][step:step+chains_per_phi_per_set],fixed_t1,'1',t2,'1')
                second_order_corr.append(soc)            
# ====            
                foc = prob_at_t(chain_all[i][step:step+chains_per_phi_per_set],t2,'1')
                first_order_corr.append(foc)
            soc_for_all_timesteps.append(second_order_corr)
            foc_for_all_timesteps.append(first_order_corr)
#====
        foc_for_all_sets.append(foc_for_all_timesteps)
        soc_for_all_sets.append(soc_for_all_timesteps)
        
#        print '*',len(chain_all[i][step:step+chains_per_phi_per_set])
        ###================= Gradient of FIRST ORDER CORR=====================
        
        # w.r.t gamma
#        gradfoc_for_all_timesteps = []        
#        for i in range(0,time_steps):
#            lst2 = [item[i] for item in foc_for_all_timesteps]
#            gradfoc_for_all_timesteps.append(np.rad2deg(np.gradient(lst2)))
#        gradfoc_for_all_sets.append(gradfoc_for_all_timesteps)
                
                # w.r.t time
        gradfoc_for_all_phis = []        
        for i in range(0,no_of_phis):
            lst2 = foc_for_all_timesteps[i]
            gradfoc_for_all_phis.append(np.rad2deg(np.gradient(lst2)))
        gradfoc_for_all_sets.append(gradfoc_for_all_phis)
        ###================= Gradient of SECOND ORDER CORR=====================
        gradsoc_for_all_timesteps = []        
        for i in range(0,time_steps):
            lst2 = [item[i] for item in soc_for_all_timesteps]
            gradsoc_for_all_timesteps.append(np.rad2deg(np.gradient(lst2)))
        gradsoc_for_all_sets.append(gradsoc_for_all_timesteps)
        
    return soc_for_all_sets, gradsoc_for_all_sets, foc_for_all_sets, gradfoc_for_all_sets

def calc_soc_grad_Adam(chain_all,total_sets,total_chains_per_phi,time_steps,no_of_phis):

    print '\r Calculating Gradient and Adams Second Order Correlations ...                                                   \r',

    chains_per_phi_per_set = int(total_chains_per_phi / total_sets)    
    fixed_t1 = 1
    soc_for_all_sets = []
    grad_for_all_sets = []
    
    for step in range(0,total_chains_per_phi,chains_per_phi_per_set): # for every sub-set
        print '\r                   Processing Set No. starting from:',step,'                                           \r',
        soc_for_all_timesteps = []
        for i in range(0,no_of_phis): # for every phi
        ###========================SECOND ORDER CORR===========================   
            print '\r Phi=',i, '      \r',      
            second_order_corr = []
            for t2 in range(0,time_steps):
                soc = sec_ord_corr_Adam(chain_all[i][step:step+chains_per_phi_per_set],fixed_t1,'1',t2,'1')
                second_order_corr.append(soc)
            soc_for_all_timesteps.append(second_order_corr)
        soc_for_all_sets.append(soc_for_all_timesteps)
        
        ###================= Gradient of SECOND ORDER CORR=====================
        grad_for_all_timesteps = []        
        for i in range(0,time_steps):
            lst2 = [item[i] for item in soc_for_all_timesteps]
            grad_for_all_timesteps.append(np.rad2deg(np.gradient(lst2)))
        grad_for_all_sets.append(grad_for_all_timesteps)
        
    return soc_for_all_sets, grad_for_all_sets    
#===============================Finding the Maximum Grad ======================
def find_max_grad(gradient_all):
    a = []
    for i in range(0,99):
        a.append(np.max(gradient_all[i]))
    timestep_of_max_grad = np.argmax(a)
    max_grad_value = np.max(a)
    phi_at_max_grad = np.argmax(gradient_all[timestep_of_max_grad])+0.5
    return (max_grad_value,phi_at_max_grad)

#===============================Plot SOC wrt Phi===============================
def plt_soc_wrt_phi(data_set,set_no,at_timestep):
    lst2 = [item[at_timestep] for item in data_set[set_no]]
    plt.plot(lst2)
    plt.axis([1, 180, -1, 2])
    
#    a = np.power( np.sin(   np.deg2rad(range(0,180))    ),4 )   
#    b = np.power( np.sin(   np.deg2rad(range(0,180))    ),2 ) 
#    c = 4 * np.sin( np.deg2rad(range(0,180) ))**3 * np.cos(np.deg2rad(range(0,180)))
#    plt.plot(a)
#    plt.plot(b)
#    plt.plot(c)    
    plt.show()

def plt_foc_wrt_phi(data_set,set_no,from_time,to_time):
        
    for t in range(from_time,to_time):
        lst2 = [item[t] for item in data_set[set_no]]
        plt.plot(lst2)
    plt.show()    
    
def plt_grad_wrt_time(data_set,set_no,time_step):
    plt.plot(data_set[set_no][time_step])
    plt.show()
#==============================================================================