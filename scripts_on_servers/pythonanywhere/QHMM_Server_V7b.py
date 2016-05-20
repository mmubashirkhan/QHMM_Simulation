# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:25:24 2016

@author: mmkhan
"""
from __future__ import division

import qutip as qt
import numpy as np
import random as rn
import multiprocessing
import pickle
import os, glob, shutil, time


class QHMMv2(multiprocessing.Process):

    def __init__(self):
        multiprocessing.Process.__init__(self)
        self.signal_value = ""
        self.state_vector = (qt.basis(2, 0) + qt.basis(2, 1)).unit() #unit method to again normalize the state. Operating with the number function again:

    def runGeneralQHMM(self, timeSteps,phi):

        self.current_state = self.state_vector
        self.signal_value = ''

# Kross Operator Adam
#        self.k0 = qt.Qobj([[np.cos(phi),0],[0,np.cos(phi)]])
#        self.k1 = qt.Qobj([[0,np.sin(phi)],[np.sin(phi),0]])

# Kross Operator Lewis
#        self.k0 = qt.Qobj([[np.cos(phi),0],[0,np.sin(phi)]])
#        self.k1 = qt.Qobj([[0,np.cos(phi)],[np.sin(phi),0]])

        del_t = 0.0001

        self.k0 = qt.Qobj([[1,0],[0, np.exp( (-1/2) * phi * del_t) ]])
        self.k1 = qt.Qobj([[0,0],[0,np.sqrt( phi * del_t ) ]])


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

def generateChains(file_name,angles):
    phi_all = []
    chain_all = []
    
    x = []
    
    for i in range(0,10):
        x.append(QHMMv2())
    
    for i in angles:
        print i
        phi_all.append(i)
    
        time_steps = 10
        no_of_chains = 100
    
        chains = []
        for j in range(0,no_of_chains):
            print j
    
            for sim in x:
                chains.append(sim.runGeneralQHMM(time_steps,np.deg2rad(i)))
    
        for c in chains:
            chain_all.append(c)
    return chain_all
#==============================================================================

#angles = np.linspace(0.5,19.5,20)  # file_a
angles = np.linspace(20.5,39.5,20)  # file_b
#angles = np.linspace(40.5,59.5,20)  # file_c
#angles = np.linspace(60.5,79.5,20)  # file_d
#angles = np.linspace(80.5,99.5,20)  # file_e
#angles = np.linspace(100.5,119.5,20)  # file_f
#angles = np.linspace(120.5,139.5,20)  # file_g
#angles = np.linspace(140.5,159.5,20)  # file_h
#angles = np.linspace(160.5,179.5,20)  # file_i
dir = str(os.getcwd()) + "/data_b_" + str(time.time())
if not os.path.exists(dir):
    os.makedirs(dir)

for i in range(0,100):
    file_name = 'file_b_v7'+ str(i)
    chain_all = generateChains(file_name,angles)
    save_list_to_file(file_name,chain_all)
    
    files = glob.iglob(os.path.join(os.getcwd(), "*.pkl"))
    for file in files:
        shutil.move(file, dir)
#==============================================================================