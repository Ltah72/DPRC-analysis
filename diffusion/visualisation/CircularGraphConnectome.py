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

#define connectome template directory
ConnectomeTemplateDir = 'V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/template/'
#go into directory
os.chdir(ConnectomeTemplateDir)

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

#define connectome data directory
ConnectomeDir = 'V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/ConnectomeWhole_linear_trend_stats/'
#go into Connectome directory
os.chdir(ConnectomeDir)

#load in the connectome matrix. 
connectome = np.loadtxt(os.path.join('outputWhole_360Nodes_connectome_linear_trend_fwe_1mpvalue_t2.csv'), delimiter=',')

#plot the circle plot. Note that the number of label names need to match the number of node colours from the atlas.  
mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0, vmax = 10, colorbar=False)
mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0.949, vmax = 1, colorbar=False)


#node_colors2[1]= [54,103,129,1]