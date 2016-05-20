# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:25:24 2016

@author: mmkhan
"""
from __future__ import division

import qutip as qt
#import scipy.linalg as spln
import numpy as np
#import math as mt
import random as rn
#import matplotlib.pyplot as plt
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
        
        self.k0 = qt.Qobj([[np.cos(phi),0],[0,np.sin(phi)]])
        self.k1 = qt.Qobj([[0,np.cos(phi)],[np.sin(phi),0]])
        
        for i in range(0, timeSteps):
            
            randomNo = rn.uniform(0, 1)    
            k0_psi = self.k0 * self.current_state
            k1_psi = self.k1 * self.current_state
    
            prob_0 = k0_psi.dag() * k0_psi            
            prob_0 = float(np.abs(prob_0[0])[0][0])            
            
            if randomNo >= np.around(prob_0, decimals=8):
                self.current_state = k0_psi.unit()
                self.signal_value += '0'
            else:
                self.current_state = k1_psi.unit()
                self.signal_value += '1'
        return self.signal_value
#===============================GENERATE DATA==================================
def generateChains(phi_start,phi_end,phi_total):

    angles = np.linspace(phi_start,phi_end,phi_total)

    phi_all = []
    chain_all = []    
    x = []
    
    for i in range(0,250):
        x.append(QHMM())
        
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
    return chain_all

#===============================Calculating probability========================
def prob_at_t(chains,t,emission_symbol):
    count = 0

    for i in range(0,len(chains)):        
        a = chains[i][t:t+len(emission_symbol)]
        b = emission_symbol

        if a == b:
            count += 1
    prob = float(count)/float(len(chains))
#    try:
#        prob = float(count)/float(len(chains))
#    except ZeroDivisionError:
#        #print 'div by zero at t = ', t
#        prob = 0.1
#        print 'foc ? '
#    print 'count:', count, 'len: ', len(chains), 'prob:', prob
    #print prob, float(count), float(len(chains))
    return prob
    
def sec_ord_corr(chains,t1,emission_symbol_1,t2,emission_symbol_2):
    #print 'new prog'
    count = 0
    foc = 0
    p = 0
    for i in range(0,len(chains)): 
        print i, t1, t2
        a_t1 = chains[i][t1]
        b_t1 = emission_symbol_1
        
        a_t2 = chains[i][t2]
        b_t2 = emission_symbol_2

        if a_t1 == b_t1:
            foc += 1
        
        if a_t2 == b_t2:
            p += 1

        if a_t1 == b_t1 and a_t2 == b_t2:
            count += 1
    print float(count), float(float(foc)) ,  float(p) , float(len(chains))
    sys.exit('testing!!!')
    soc = (float(count)/float(float(foc))) /  (float(p) / float(len(chains)))

    return soc
#============================saving data to a file=============================
def save_list_to_file(file_name,list_name):
    with open(file_name+'.pkl', 'wb') as pickle_file:
        pickle.dump(list_name, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
    return 0
#============================Extract Data======================================        
def restore_list_from_file(filename_list):
    new_list = []
    for file_name in filename_list:
        print file_name
        with open('Sim_data_final/'+file_name+'.pkl', 'rb') as pickle_load:
            temp_list = pickle.load(pickle_load)
        new_list.append(temp_list)
    return new_list  
#========================Reformatting/Gathering Data===========================
def extract_chains_per_phi(lst,chains_per_angle):
    lst_a = []
    a = 0
    d = chains_per_angle
    for n in range(1,361):
        term = a + (n-1)*d         
        lst_a.append(lst[term:term+chains_per_angle])
    return lst_a
#==============================================================================