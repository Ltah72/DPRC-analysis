# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:22:55 2021

@author: ltah262

%Author: Lenore Tahara-Eckl
%Email: Ltah262@aucklanduni.ac.nz
%Date: 1/02/21

Use the python package, MNE to create a circle plot for the connectome data. 
https://mne.tools/stable/index.html
"""

#Import necessary libraries
import os
import numpy as np
import mne

#define connectome data directory
ConnectomeDir = 'V:/NECTAR_data/LENORE/test/connectome_test/'
#go into Connectome directory
os.chdir(ConnectomeDir)

#Match atlas colours to the node names
labels = []
node_colors = []
unknown_found = None
with open('hcpmmp1_ordered_360nodes.txt') as f:
        for x in f.readlines():
            if unknown_found and x.split():
                labels.append(x.split()[1]) 
                node_colors.append(tuple([float(color)/255 for color in x.split()[2:6]]))
            if 'Unknown' in x:
                unknown_found = True

#Create the circular layout by determining the number of angles in the circle, based upon the number of node names. 
#Also, divide and group the circle plot into the left and right hemisphere. 
node_angles = mne.viz.circular_layout(labels, labels, group_boundaries=[0, len(labels) / 2])

#load in the connectome matrix. 
connectome = np.loadtxt(os.path.join('outputWhole_connectome_fwe_1mpvalue_t4_360nodes.csv'), delimiter=',')

#plot the circle plot. Note that the number of label names need to match the number of node colours from the atlas.  
#mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors)
mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0.94, vmax = 1, colorbar=False)
