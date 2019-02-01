#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 10:32:40 2018
Script to read a mat file and analyze the particle tracking script 
@author: apreziosir
"""

import numpy as np
import h5py
import time

def Conc_map(nameI, nameO):
    
    # Timing the operation (just for exploration, it can be optimized)
    # t0 = time.time()
    # print('Empiezo a las: \t', t0)

    # ==========================================================================
    # Defining variables for the script (mainly dx and dy)
    # ==========================================================================
    
    # Cell size in non dimensional form 
    dx = np.pi / 20
    dy = np.pi / 20
    
    # ==========================================================================
    # HANDLING THE HDF5 FILE, EXTRACTING VALUES AND RUNNING IT
    # ==========================================================================
    
    # Read an hdf5 file and store it in the x variable. The file must be read 
    # and the values extracted before closing the hdf5 file that stored the 
    # results from the model. Otherwise is not useful
    f = h5py.File(nameI, 'r')
    x = f['DATA_1']
    
    # Loading the complete positions (moving andd filtered) - loading the values
    # just once, so it is not necessary to access the hdf5 every timestep
    x_mov = x['keepx'][:,:]
    y_mov = x['keepy'][:,:]
    x_fil = x['keepxfil'][:,:]
    y_fil = x['keepyfil'][:,:]
    
    # Extracting values of the hdf5 file (dimension 1 x 1)
    Lx = x['Lx'][0]
    Ly = x['Ly'][0]
    n = x['n'][0]
    
    # Extracting 1D arrays from the HDF5 file
    keepframes = x['keepframes'][:]
    
    # Closing the hdf5 file from which I read data
    f.close()
    
    # ==========================================================================
    # Defining variables that depend on program inputs and data read from the 
    # hdf5 file 
    # ==========================================================================
    
    # Numer of cells in the domain that is considered, not floor, not ceil
    Nx = np.round(Lx / dx, decimals = 0).astype(int)
    Ny = np.round(Ly / dy, decimals = 0).astype(int)
    
    # Coordinates of the mesh in vector form
    cX = np.linspace(0, Lx[0], Nx[0] + 1)
    cY = np.linspace(-Ly[0], 0, Ny[0] + 1)
    
    # Creating 3D arrays that store the information of the map in "layers"
    conc_map = np.empty((Ny[0], Nx[0]))
    filt_map = np.empty((Ny[0], Nx[0]))
    tot_conc = np.empty((Ny[0], Nx[0]))
    
    # ==========================================================================
    # Creating an hdf5 file to store the results of the script. Groups are also 
    # created in order to store things in an efficient and readable way
    # ==========================================================================
    
    # Creatng file - We have to create file in 'w' mode and then close it and 
    # re-open it in 'r+' mode
    f1 = h5py.File(nameO, 'w')
    f1.close()
    
    WRT = h5py.File(nameO, 'a')
    
    # Creating file groups to store in different variables
    Gmvng = WRT.create_group('Moving_part')
    Gfilt = WRT.create_group('Filter_part')
    Gtotal = WRT.create_group('Total_part')
    
    # ==========================================================================
    # Iterating over time and space to get the relative concentration maps
    # First row of filtered positions is filled with 0, so the loop must start 
    # at index 1 (second row)
    # ==========================================================================
    
    # Variable to check if all particles have stopped or have been remobilized 
    check_x = 0
    check_y = 0
    
    for t in range(1, len(keepframes)):
        
        # Timing the process (just to check how log it takes to search 
        # particles)
    #    t1 = time.time()
    #    print('Elapsed time in iteration: ', t, '\t', t1 - t0)
        
        # Extracting matrices with particles' position for every timestep saved
        # This operation is performed once for every timestep so it is quicker 
        # than setting it up inside the three loops
        posx = x_mov[t, :]
        posy = y_mov[t, :]
        
        posxfil = x_fil[t, :]
        posyfil = y_fil[t, :]   
        
        # Names of datasets for saving hdf5 file
        c_name = 'Conc_' + str(t).zfill(4)
        f_name = 'Filt_' + str(t).zfill(4)
        t_name = 'Totl_' + str(t).zfill(4)
        
        for j in range(0, len(cY) - 1): 
           
            for i in range(0, len(cX) - 1):
                    
                # Looking for particles in each cell. Since the keepx and keepy
                # arrays have the same size, we can mix conditions (filtered and 
                # moving particles)
                parts = np.where(np.logical_and(np.logical_and(posx >= cX[i], \
                        posx < cX[i + 1]), np.logical_and(posy >= cY[j], posy \
                        < cY[j + 1])))
                
                partsf = np.where(np.logical_and(np.logical_and(posxfil >= cX[i], \
                     posxfil < cX[i + 1]), np.logical_and(posyfil >= cY[j], \
                     posyfil < cY[j + 1])))
                
                # In this place the code shows a warning using comparison 
                # operators like <, <=, >, >=. This is due to the presence of 
                # nan in the input arrays. However, it seems to work OK
                
                # Making the vectors found numpy arrays
                parts = np.asarray(parts)
                partsf = np.asarray(partsf)
                
                # Printing size of vectors for checking
    #            print('Size of arrays: ')
    #            print(parts.size, partsf.size)
                
                # The vector lengths of part and partsf are the number of 
                # particles in each cell
                conc_map[j, i] = parts.size
                filt_map[j, i] = partsf.size
               
        # Adding the concentration values in the whole three dimensional array 
        # (this is outside the temporal loop since filtered and moving are in 
        # separate loops)
        tot_conc = conc_map + filt_map
        
        # Saving the 2D array into the hdf5 file to be manipulated later
        dset1 = Gmvng.create_dataset(c_name, data = conc_map)
        dset2 = Gfilt.create_dataset(f_name, data = filt_map)
        dset3 = Gtotal.create_dataset(t_name, data = tot_conc)
        
        # Checking if all the particles have been filtered or remobilized
        check_x = np.sum(np.isnan(posx))
        check_y = np.sum(np.isnan(posy))
        
        if check_x == n and check_y == n:
            
            break
                
    # Closing the hdf5 files. This operation is done when all the values from 
    # the hdf5 file are extracted into variables. Otherwise they won't be saved.
    # THIS IS VERY VERY VERY VERY IMPORTANT!!!! EXTREMELY IMPORTANT!!!!
    WRT.close()
    
    # Timing (exploration purposes)
    # t1 = time.time() - t0
    # print('termino a las: \t', t1)
    # print('Elapsed time: ', t1, ' seconds')
    
    return()
