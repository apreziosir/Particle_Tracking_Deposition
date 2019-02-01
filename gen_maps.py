#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 08:50:35 2018
Script to create file names and generate hdf5 files of the particle tracking
script. 
@author: Antonio Preziosi-Ribero
"""

import os
import readmatF as RF

# Creating vector of in filenames, which will work also as generator of 
# out names for the hdf5 files
in_name0 = '190108_FinBed/'
out_folder = '190116_HDF5/'
IN = os.listdir(in_name0)

# Counting the number of times the process has been executed
count = 0

# Iterating over arrays, generate filename out
for II in IN:
    
    # In file name (since it is in other folder, I have to change it)
    in_name = in_name0 + II
#    print(in_name)             # Checking name consistency
    
    # Generate name of the HDF5 file that is going to be saved
    out_name = out_folder + 'Map_' + II[5:-4] + '.hdf5'
#    print(out_name)            # Checking name consistency
    
    # Running function that will create hdf5 file
    RF.Conc_map(in_name, out_name)
    
    # Printing log in screen
    print('Just finished iteration number: \t', count + 1)
    print('Out of: \t', len(IN))
    print('Advance: \t', ((count + 1) / len(IN)) * 100, '%')
    
    # Advancing counter
    count += 1   
