#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 16:58:19 2025

For each simulation use the respective parameter:
    
 GM (Control case A): parameter = 'A'  
 GZ+GM (Control case B): parameter = 'B'
 PUSH
    gamma - parameter
    0.05 - '5em2'
    0.1 - '10em2'
    0.2 - '20em2'
    0.3 - '30em2'
    0.5 - '50em2'
    1.0 - '1'
    2.0 - '2'
    3.0 - '3'
 PULL   
    gamma_hat - parameter
    0.05 - '5em2'
    0.1 - '10em2'
    0.2 - '20em2'
    0.3 - '30em2'
    0.5 - '50em2'
    1.0 - '1'
    2.0 - '2'
    3.0 - '3'

This generates .csv files that are read by the main code to get the overall GI trends.
This also generates accompanying individual trends for GI for time. 
@author: ktaneja
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.collections import LineCollection
from Gyrification_functions import *
from pathlib import Path
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'#'/latex'

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 50
LARGE_SIZE = 18

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

plt.rcParams["font.family"] = "Times New Roman"

current_folder = Path(__file__).resolve().parent ## 'sim-analysis-GI'
main_directory = current_folder.parent ## 'simulation-files'
odb_data_folder = main_directory/"odb-data-files"


model_name_flag = 'control' #'PULL', 'control' 
parameter = 'B' #'30em2' #'5em2','10em2','20em2','50em2','A','B' 


if model_name_flag=='control':

    cortex_top_coord_disp_file = odb_data_folder/"cortex_top_coords_disp_control{}.csv".format(str(parameter))
    
    hull_fig_filename = 'convex_hull_file_control{}.pdf'.format(str(parameter))
    fig_filename = 'GI_file_control{}.pdf'.format(str(parameter))
    result_filename = 'GI_control{}.csv'.format(str(parameter))

    
if model_name_flag=='PUSH':

    cortex_top_coord_disp_file = odb_data_folder/"cortex_top_coords_disp_gamma{}.csv".format(str(parameter))
     
    hull_fig_filename = 'convex_hull_gamma{}.pdf'.format(str(parameter))
    fig_filename = 'GI_gamma{}.pdf'.format(str(parameter))
    result_filename = 'GI_gamma{}.csv'.format(str(parameter))


if model_name_flag=='PULL':

    cortex_top_coord_disp_file = odb_data_folder/"cortex_top_coords_disp_gammahat{}.csv".format(str(parameter))
    hull_fig_filename = 'convex_hull_file_gammahat{}.pdf'.format(str(parameter))
    fig_filename = 'GI_file_gammahat{}.pdf'.format(str(parameter))
    result_filename = 'GI_gammahat{}.csv'.format(str(parameter))



save_dir_figures = "figures"
# save_dir_GI = "GI_sims"
os.makedirs(save_dir_figures, exist_ok=True)  # Create directory if it doesn't exist
# os.makedirs(save_dir_GI, exist_ok=True)  # Create directory if it doesn't exist


time_range = [0.1,0.15,0.2
              ,0.25,0.3,0.35,
              0.4,0.42,0.45,0.47,
              0.5,0.51,0.52,0.55,
              0.6,0.62,0.65,
              0.7,0.72,0.75,0.77,
              0.8,0.82,0.85,0.87,
              0.9,0.92,0.97,0.99]

 
#%%
        
## GM processing
## Import coords and dispm file
cortex_top_coord_disp_df = pd.read_csv(cortex_top_coord_disp_file)

z_slice = 0.6#1.0#0.6669773
xmin = 0.1 # Consider points slightly off the boundary
ymin = 0.1 # Consider points slightly off the boundary

# Consider points slightly off the boundary and on a slice at cordiante z_slice. 
tol = 1e-3
cortex_top_coord_disp_df = cortex_top_coord_disp_df[cortex_top_coord_disp_df['X']>xmin]
cortex_top_coord_disp_df = cortex_top_coord_disp_df[cortex_top_coord_disp_df['Y']>ymin]
cortex_coord_disp_zslice_df = cortex_top_coord_disp_df[abs(cortex_top_coord_disp_df['Z']-z_slice)<tol]



GM_columns = ['Node','X','Y','Z','U1','U2','U3']
cortex_top_disp_dict = {time: None for time in time_range}


## Get deformed coordinates of cortex
for time in time_range:
    num_rows = cortex_coord_disp_zslice_df.shape[0]   
    df = pd.DataFrame(columns=GM_columns, index=range(num_rows))
    for label in range(cortex_coord_disp_zslice_df.shape[0]):
        node_ID = cortex_coord_disp_zslice_df.iloc[label,0]
        
        ref_coords2 = cortex_coord_disp_zslice_df[cortex_coord_disp_zslice_df['Node']==node_ID].iloc[:,1:4].to_numpy()
    ## Get the columns for current time
        node_mask = cortex_coord_disp_zslice_df.columns.str.endswith(str(time))
        deformations = cortex_coord_disp_zslice_df[cortex_coord_disp_zslice_df['Node']==node_ID].loc[:,node_mask].to_numpy()
    ## Store in df
        df.iloc[label,0] = node_ID
        df.iloc[label,1:4] = ref_coords2
        df.iloc[label,4:] = deformations
    ### Store df in dict
        cortex_top_disp_dict[time] = df
       
#%%
               

def sort_contour_points(points):
    # Compute centroid
    center = np.mean(points, axis=0)

    # Compute angles of each point w.r.t. the centroid
    angles = np.arctan2(points[:, 1] - center[1], points[:, 0] - center[0])

    # Sort points by angle
    sorted_indices = np.argsort(angles)
    return points[sorted_indices],sorted_indices
    
#%%
                  
disp_gray_data_list = []


for t in time_range:
    disp_gray_data_list.append(cortex_top_disp_dict[t])
#%%
L_oc_list = []
L_ch_list = []
GI_list = []


    
#%%

figsize = (8,2*len(disp_gray_data_list))


ms = 5-1
lw = 5-1
widths = np.array([0.005])
fig,axs = plt.subplots(len(disp_gray_data_list),figsize=figsize)

for file_ID in range(len(disp_gray_data_list)):
    print(file_ID)
    df_gray = disp_gray_data_list[file_ID]

    # Get points for convex Hull
    ref_coords_top = []
    def_coords_top = []
    ref_coords_bottom = []
    def_coords_bottom = []

## Get coordiantes of pial surface and gray-white matter interface
    for j in range(df_gray.shape[0]):
        data = df_gray.iloc[j].to_numpy()
        ref_coord = [data[1], data[2]]
        ref_coords_top.append([ref_coord[0], ref_coord[1]])
        def_coords_top.append([ref_coord[0] + data[4], ref_coord[1] + data[5]])
      

    def_coords_top = np.array(def_coords_top)
    ref_coords_top = np.array(ref_coords_top)

    top_coords_sorted_ref,top_indices_sorted_ref = sort_contour_points(ref_coords_top)   

    gray_surface_points = def_coords_top[top_indices_sorted_ref]        
       
    axs[file_ID].plot(gray_surface_points[:,0],gray_surface_points[:,1],
             '-',lw=lw,color='lightgrey',label='Pial surface')
  
    gray_surface_hull = ConvexHull(gray_surface_points)
    ## Some points are below the convex hull , so need to remove them.
    gray_surface_hull_vertices = gray_surface_hull.vertices

            
    gray_surface_hull_vertices = np.sort(gray_surface_hull_vertices)
    gray_surface_hull_coords = gray_surface_points[gray_surface_hull_vertices]
    x = gray_surface_hull_coords[:, 0]
    keep = np.ones(len(x), dtype=bool)
    keep[1:-1] &= ~((x[1:-1] < x[:-2]) & (x[1:-1] < x[2:]))
    gray_surface_hull_coords = gray_surface_hull_coords[keep]   
    
    axs[file_ID].plot(gray_surface_hull_coords[:,0],
             gray_surface_hull_coords[:,1], 'r--',markersize=ms, lw=lw,label='Convex hull')
    axs[file_ID].set_yticks([])
    axs[file_ID].set_xticks([])
    axs[file_ID].axis('off')
    axs[file_ID].set_aspect('equal')

    # Positioning legends below the plot
    plt.legend(loc='upper center',
               bbox_to_anchor=(0.5, 0.0),
               ncol=1,
               frameon=True)
    # Calculate GI

    L_oc = polyline_length(gray_surface_points)
    L_ch = polyline_length(gray_surface_hull_coords)    

    
    GI = L_oc/L_ch
    GI_list.append(GI)
    L_oc_list.append(L_oc)
    L_ch_list.append(L_ch)
    print('GI=',GI)


#%%
lw = 3
fig, ax = plt.subplots()
# --- guide lines at min y and min x from gray_surface_points ---
# horizontal line
y_min_idx = np.argmin(gray_surface_points[:, 1])
y_min = gray_surface_points[y_min_idx, 1]
x_max = gray_surface_points[y_min_idx, 0]

x_min_idx = np.argmin(gray_surface_points[:, 0])
x_min = gray_surface_points[x_min_idx, 0]
y_max = gray_surface_points[x_min_idx, 1]

# horizontal line at lowest y
# ax.axhline(y_min, linewidth=lw, color='k', zorder=5)
ax.plot([x_min,x_max],[y_min,y_min], linewidth=lw, color='k',alpha=0.7)

# vertical line at lowest x
# ax.axvline(x_min, linewidth=lw, color='k', zorder=5)
ax.plot([x_min,x_min],[y_min,y_max], linewidth=lw, color='k',alpha=0.7)

ax.plot(gray_surface_points[:, 0], gray_surface_points[:, 1], '-', lw=lw, color='lightgrey', label='Pial surface')
ax.plot(gray_surface_hull_coords[:, 0], gray_surface_hull_coords[:, 1],
        'r--', lw=lw, label='Convex hull')  





ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal', adjustable='box')
ax.axis('off')
plt.savefig(os.path.join(save_dir_figures, hull_fig_filename),format='pdf')

#%%
# Create a single-column subplot layout
# fig, axs = plt.subplots(2, 1, figsize=(20, 80))
fig, axs = plt.subplots(1, 1, figsize=(10, 20))

# GI plot
axs.plot(time_range, GI_list, '-*', linewidth=lw)
axs.set_ylabel("Gyrification Index (GI)")
axs.set_xlabel("Simulation Time")
axs.grid(True)
axs.set_ylim(1.0, 1.6)


plt.tight_layout()

# Save the figure

plt.savefig(os.path.join(save_dir_figures, fig_filename),format='pdf')


AAcopy_data = np.column_stack((GI_list,time_range))

header = "GI,SimTime"
# np.savetxt(os.path.join(save_dir_GI, result_filename), AAcopy_data, delimiter=',', header=header)
np.savetxt( result_filename, AAcopy_data, delimiter=',', header=header)
