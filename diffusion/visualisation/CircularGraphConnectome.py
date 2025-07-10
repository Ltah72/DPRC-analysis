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
ConnectomeTemplateDir = 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/template/'
#go into directory
os.chdir(ConnectomeTemplateDir)

#Match atlas colours to the node names
labels = []
node_colors = []
unknown_found = None
with open('fs_default_ordered_circular_plot_input.txt') as f:
        for x in f.readlines():
            if unknown_found and x.split():
                labels.append(x.split()[1]) 
                node_colors.append(tuple([float(color)/255 for color in x.split()[2:6]]))
            if 'Unknown' in x:
                unknown_found = True

#Create the circular layout by determining the number of angles in the circle, based upon the number of node names. 
#Also, divide and group the circle plot into the left and right hemisphere. 
#node_angles = mne.viz.circular_layout(labels, labels, group_boundaries=[0, len(labels) / 2])
#Divide the circular plot layout by cortical, cerebellum, and subcortical regions
#node_angles = mne.viz.circular_layout(labels, labels, group_boundaries = [0, 34, 68, 70, 77])
node_angles = mne.viz.circular_layout(labels, labels, group_boundaries = [0, 34, 35, 42, 76, 77])


#define connectome data directory
ConnectomeDir = 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/thresholded/wgt_strmln_cns_thr/ConnectomeWhole_stats_dk_linear_trend_noAD_cov-age_sex/'
#ConnectomeDir = 'V:/Archive/NECTAR_data/LENORE/derivatives/groups/F2/diff_data/longitudinal/connectome/stats_results/weighted/thresholded/C_SCD_threshold/ConnectomeWhole_stats_dk_linear_trend_cov-age_sex/'

#go into Connectome directory
os.chdir(ConnectomeDir)

#load in the connectome matrix. 
#connectome = np.loadtxt(os.path.join('transferred_values_big-node_linear_trend_82_nodes_for_circle_plot_FPN_stats-fwe_1mpvalue_t1.csv'), delimiter=',')
connectome = np.loadtxt(os.path.join('pval_display_CircPlot_linear_trend_noAD_fwe_1mpvalue_t1.csv'), delimiter=',')

#plot the circle plot. Note that the number of label names need to match the number of node colours from the atlas.  
#mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0, vmax = 10, node_edgecolor='white', colorbar=False)
mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0.94, vmax = 1, colormap = 'hot', colorbar=False)

#mne.viz.plot_connectivity_circle(connectome, labels, node_angles = node_angles, node_colors = node_colors, vmin = 0.95, vmax = 0.96, colormap = 'hot', colorbar=False)


#node_colors2[1]= [54,103,129,1]