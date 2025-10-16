#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:32:10 2024

@author: ktaneja
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = "Arial"

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

 
dpi = 600*2*2 # resolution for figures

## For stress figures
colorbar_name = 'Colorbar_NegOne_to_PosOne.pdf'

fig, ax = plt.subplots(figsize=(1, 8)) ## Vertical orientation
norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0) ## Fig 9,11: Stress colorbar
vals = [-1.0,-0.5,0.0,0.5,1.0] ## Fig 9,11: Stress colorbar ticks
cmap = plt.cm.coolwarm_r ## Stress

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')


## For every other figure
# colorbar_name = 'Colorbar_One_to_Three.pdf'
# fig, ax = plt.subplots(figsize=(8, 1)) ## Horizontal orientation
# cmap = mpl.cm.plasma ## Everything else
# norm = mpl.colors.Normalize(vmin=1.0, vmax=3.0) ## Fig5-7, 10: Growth colorbar
# vals = [1.0,1.5,2.0,2.5,3.0] ## Fig5-7, 10: Growth colorbar ticks



# norm = mpl.colors.Normalize(vmin=0.0, vmax=1.3) ## Fig 2 A colorbar
# vals = [0.0,0.3,0.65,1.0,1.3] ## Fig 2 A colorbar ticks

# cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
#                                 norm=norm,
#                                 orientation='horizontal')

cb1.set_ticks(vals,labels=vals)

plt.tight_layout()

plt.savefig(colorbar_name,format='pdf',dpi=dpi,bbox_inches='tight')