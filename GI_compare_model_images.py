#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 22:24:11 2025

@author: ktaneja
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import splrep
from scipy.interpolate import splev
import seaborn as sns
## Sulcal Depth
import math
from scipy.spatial import ConvexHull

import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'#'/latex'

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 100#40+5+5+5+5+5+5
LARGE_SIZE = 18

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

# plt.rcParams["font.family"] = "Times New Roman"
# Set Palatino as the default font
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Palatino"]
plt.rcParams['xtick.major.pad'] = 15
plt.rcParams['ytick.major.pad'] = 20
plt.rcParams['xtick.minor.pad'] = 12
plt.rcParams['ytick.minor.pad'] = 12
# Enable LaTeX rendering
# plt.rc('text', usetex=True)
# plt.rcParams['text.usetex'] = True
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.serif'] = ['Computer Modern Roman'] # Or other desired font
# plt.rcParams['text.latex.preamble'] = r'\usepackage{lmodern}' # Example for Latin Modern


ms = 15+5+5 # markersize
lw = 7+2 # linewisdth
color_thresh = -0.1
save_dir_figures = "manuscript-figures"
os.makedirs(save_dir_figures, exist_ok=True)  # Create directory if it doesn't exist
legend_flag = 1

axlim = (1.0,16.0)
lp = 20+5 # Label padding from axes
fsize = (20, 25+5)
alpha_level = 0.4-0.1  # Adjust transparency (0: fully transparent, 1: opaque)

ticks = [1,4,8,12,16]
color_thresh = -0.1

P10_GI = [1.035,1.093]
Age_10 = [10,10]

Age_16 = [16,16,16]  
P16_GI = [1.401072817,1.340624800701704,1.2536]


# --- helper to load a mode (unchanged idea) ---
def load_mode(mode):
    assert mode in ('push','pull')

    if mode == 'pull':
        result_file_gammahat0p05  = 'simulation-files/sim-analysis-GI/GI_gammahat5em2.csv'
        result_file_gammahat0p1 = 'simulation-files/sim-analysis-GI/GI_gammahat10em2.csv'
        result_file_gammahat0p2 = 'simulation-files/sim-analysis-GI/GI_gammahat20em2.csv'
        result_file_gammahat0p3 = 'simulation-files/sim-analysis-GI/GI_gammahat30em2.csv'
        result_file_gammahat0p5 = 'simulation-files/sim-analysis-GI/GI_gammahat50em2.csv'
        result_file_gammahat1= 'simulation-files/sim-analysis-GI/GI_gammahat1.csv'
        result_file_gammahat2= 'simulation-files/sim-analysis-GI/GI_gammahat2.csv'
        result_file_gammahat3= 'simulation-files/sim-analysis-GI/GI_gammahat3.csv'
        result_file_control_osvz= 'simulation-files/sim-analysis-GI/GI_controlB.csv'
        result_file_control_crtx= 'simulation-files/sim-analysis-GI/GI_controlA.csv'

        AA = [
            np.genfromtxt(p, delimiter=',', skip_header=1, usecols=(0,1)).reshape((-1,2))
            for p in [
                result_file_control_crtx, result_file_control_osvz,
                result_file_gammahat0p05, result_file_gammahat0p1, result_file_gammahat0p2,
                result_file_gammahat0p3, result_file_gammahat0p5, result_file_gammahat1,
                result_file_gammahat2, result_file_gammahat3
            ]
        ]
        GI, SimTime = [[a[:,i] for a in AA] for i in range(2)]
        labels = [
            'Control (GM)', r'Control (GZ + GM)',
            r'$\gamma,\hat{\gamma}=0.05$', r'$\gamma,\hat{\gamma}=0.1$', r'$\gamma,\hat{\gamma}=0.2$', r'$\gamma,\hat{\gamma}=0.3$',
            r'$\gamma,\hat{\gamma}=0.5$',  r'$\gamma,\hat{\gamma}=1.0$', r'$\gamma,\hat{\gamma}=2.0$', r'$\gamma,\hat{\gamma}=3.0$'
        ]
        t = np.linspace(1, 16, num=AA[0].shape[0])
        return dict(GI=GI, SimTime=SimTime, labels=labels, t=t)

    else:  # push
        result_file_gamma0p05 = 'simulation-files/sim-analysis-GI/GI_gamma5em2.csv'
        result_file_gamma0p1  = 'simulation-files/sim-analysis-GI/GI_gamma10em2.csv'
        result_file_gamma0p2  = 'simulation-files/sim-analysis-GI/GI_gamma20em2.csv'
        result_file_gamma0p3  = 'simulation-files/sim-analysis-GI/GI_gamma30em2.csv'
        result_file_gamma0p5  = 'simulation-files/sim-analysis-GI/GI_gamma50em2.csv'
        result_file_gamma1    = 'simulation-files/sim-analysis-GI/GI_gamma1.csv'
        result_file_gamma2    = 'simulation-files/sim-analysis-GI/GI_gamma2.csv'
        result_file_gamma3    = 'simulation-files/sim-analysis-GI/GI_gamma3.csv'
        result_file_control_osvz= 'simulation-files/sim-analysis-GI/GI_controlB.csv'
        result_file_control_crtx= 'simulation-files/sim-analysis-GI/GI_controlA.csv'

        AA = [
            np.genfromtxt(p, delimiter=',', skip_header=1, usecols=(0,1)).reshape((-1,2))
            for p in [
                result_file_control_crtx, result_file_control_osvz,
                result_file_gamma0p05, result_file_gamma0p1, result_file_gamma0p2,
                result_file_gamma0p3, result_file_gamma0p5, result_file_gamma1,
                result_file_gamma2, result_file_gamma3
            ]
        ]
        GI, SimTime = [[a[:,i] for a in AA] for i in range(2)]
        labels = [
            'Control (GM)', r'Control (GZ + GM)',
            r'$\gamma=0.05$', r'$\gamma=0.1$', r'$\gamma=0.2$', r'$\gamma=0.3$',
            r'$\gamma=0.5$',  r'$\gamma=1.0$', r'$\gamma=2.0$', r'$\gamma=3.0$'
        ]
        t = np.linspace(1, 16, num=AA[0].shape[0])
        return dict(GI=GI,SimTime=SimTime, labels=labels, t=t)


# --- build 1x2 grid: GI only (push | pull) ---
# save_dir_figures = "figures"
# os.makedirs(save_dir_figures, exist_ok=True)

push = load_mode('push')
pull = load_mode('pull')

fig, (ax_gi_push, ax_gi_pull) = plt.subplots(
    1, 2, figsize=(fsize[0]*3.0, fsize[1]*0.8), constrained_layout=True
)

cmap_mode = {'push': 'plasma_r', 'pull': 'plasma_r'}

def plot_gi(ax_gi, data, mode_name):
    cmap = cm.get_cmap(cmap_mode[mode_name], len(data['GI']))
    t = data['t']
    handles = []
    labels  = []
    for k, arr in enumerate(data['GI']):
        color = cmap(k / max(1, len(data['GI'])-1) - color_thresh)
        h, = ax_gi.plot(t, arr, linewidth=lw, color=color, label=f"{data['labels'][k]}")
        handles.append(h)
        labels.append(f"{data['labels'][k]}")
    return handles, labels

# plot GI for push and pull
handles_push, labels_push= plot_gi(ax_gi_push, push, 'push')
handles_pull, labels_pull = plot_gi(ax_gi_pull, pull, 'pull')

# overlay P10/P16 markers on both GI panels
for ax in (ax_gi_push, ax_gi_pull):
    ax.plot(Age_10, P10_GI, 's', ms=25, color='dodgerblue')
    ax.plot(Age_16, P16_GI, 's', ms=25, color='dodgerblue')

# cosmetics
for ax in (ax_gi_push, ax_gi_pull):
    ax.set_xlabel(r"Postnatal age [day]")
    ax.set_xticks(ticks)
    ax.grid(True, axis='y')

ax_gi_push.set_ylabel(r"GI [-]", labelpad=lp)
ax_gi_pull.set_ylabel(r"GI [-]", labelpad=lp)

# y-limits shared for both
ax_gi_push.set_ylim(1.0, 1.6)
ax_gi_pull.set_ylim(1.0, 1.6)

# optional column titles
# ax_gi_push.set_title('Push')
# ax_gi_pull.set_title('Pull')

# single, shared legend built from the PULL GI handles only
if legend_flag:
    fig.legend(
        handles=handles_pull,
        labels=labels_pull,
        loc='center left',
        bbox_to_anchor=(0.95, 0.5),
        ncol=1,
        frameon=False
    )

plt.tight_layout()
plt.subplots_adjust(
    left=0.2, right=0.9, bottom=0.12, top=0.9, wspace=0.25, hspace=0.0
)

fig_filename = 'Fig8_GI_push_pull_ferret.pdf'
plt.savefig(os.path.join(save_dir_figures, fig_filename), format='pdf', bbox_inches='tight')
