# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:25:24 2016

@author: mmkhan
"""
import os, glob, shutil, time

def moveFiles():
    files = glob.iglob(os.path.join(os.getcwd(), "*.pkl"))
    dir = str(os.getcwd()) + "/data" + str(time.time())
    if not os.path.exists(dir):
        os.makedirs(dir)
    for file in files:
        shutil.move(file, dir)
moveFiles()