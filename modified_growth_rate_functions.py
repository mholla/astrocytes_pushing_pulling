#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 07:57:00 2024

@author: ktaneja
To be used in the figure for papers
For Latex in Mac, https://stackoverflow.com/questions/70594392/failed-to-process-string-with-tex-and-i-have-latex
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
import matplotlib as mpl
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'#'/latex'



SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 50+10+10 # 150 for original size
LARGE_SIZE = 24
lp = 10 # label padding
flag='phi' # H  else  phi
# flag='H' # H  else  phi
fig_filename = 'Fig1_phase3_growth_rate_{}_modified.pdf'.format(flag)
dpi = 600*2*2 # resolution for figures


plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Enable LaTeX rendering
# plt.rc('text', usetex=True)


# plt.rcParams["font.family"] = "Times New Roman"
Cortical_Thickness = 1.5 # NOTE: For figure to show a visible thickness, can change this to 3. (1/10 mm)
MajorAxis_G = 36  # Grey major axis in 1/10 mm
MinorAxis_G = 30  # Grey minor axis in 1/10 mm
MajorAxis_W = MajorAxis_G - Cortical_Thickness # White major axis
MinorAxis_W = MinorAxis_G - Cortical_Thickness  # White minor axis

threshold = 10
minor_reduced = 4
maj_min_ratio = MajorAxis_W/MinorAxis_W
major_reduced = minor_reduced*maj_min_ratio
MajorAxis_T = MajorAxis_W - threshold
MinorAxis_T = MinorAxis_W - threshold
thresh = 0.6 # Trying parameter to approximate VZ and germinal zone thickness.
phi_crit = np.pi / 14
phi_gap = np.pi / 14
N_gyri = 4
lw = 10 #linewidth

# =============================================================================
# USER-DEFINED SMOOTHING PARAMETERS
# =============================================================================
# SMOOTHING_TYPE: 'sigmoid' or 'smoothstep' or 'tanh'
SMOOTHING_TYPE = 'sigmoid'  # Choose: 'sigmoid', 'smoothstep', 'tanh'

# SMOOTHING_STEEPNESS: Controls how sharp the transition is
# Higher values = sharper transition (closer to step function)
# Lower values = smoother transition
SMOOTHING_STEEPNESS = 20  # Typical range: 10-50

# TRANSITION_WIDTH: For smoothstep method only
# Width of the transition region (fraction of r)
TRANSITION_WIDTH = 0.15  # Typical range: 0.05-0.3
# =============================================================================

def y_coord(x_coord, majoraxis, minoraxis):
    output = np.sqrt((1 - (x_coord * x_coord) / (majoraxis * majoraxis))) * minoraxis
    return output

def x_coord(y_coord, majoraxis, minoraxis):
    output = np.sqrt((1 - (y_coord * y_coord) / (minoraxis * minoraxis))) * majoraxis
    return output

def heaviside(x, threshold, gamma=10):
    delta = threshold
    output = np.exp(gamma * ((x - delta))) / (1 + np.exp(gamma * ((x - delta))))
    return output

def gaussian(x, mean, std):
    fac = 1 / (std * np.sqrt(2 * np.pi))
    scaled_input = (x - mean) / std
    output = fac * np.exp(-0.5 * (scaled_input ** 2))
    return output

def smoothstep_transition(r, threshold, width):
    """Smoothstep function for C1 continuous transition"""
    t = (r - (threshold - width/2)) / width
    t = np.clip(t, 0, 1)
    return t**3 * (t * (t * 6 - 15) + 10)

def sigmoid_transition(r, threshold, steepness):
    """Sigmoid (logistic) transition"""
    return 1 / (1 + np.exp(-steepness * (r - threshold)))

def tanh_transition(r, threshold, steepness):
    """Hyperbolic tangent transition"""
    return 0.5 * (1 + np.tanh(steepness * (r - threshold)))

def apply_smooth_transition(r, threshold, value_below, value_above):
    """
    Apply smooth transition between value_below (r < threshold) and 
    value_above (r > threshold) using the user-selected method
    """
    if SMOOTHING_TYPE == 'sigmoid':
        transition = sigmoid_transition(r, threshold, SMOOTHING_STEEPNESS)
    elif SMOOTHING_TYPE == 'tanh':
        transition = tanh_transition(r, threshold, SMOOTHING_STEEPNESS)
    elif SMOOTHING_TYPE == 'smoothstep':
        transition = smoothstep_transition(r, threshold, TRANSITION_WIDTH)
    else:
        raise ValueError(f"Unknown smoothing type: {SMOOTHING_TYPE}")
    
    return value_below * (1 - transition) + value_above * transition

def scaled_xy(x,y):
    MajorAxisLength = MajorAxis_W - major_reduced
    MinorAxisLength = MinorAxis_W - minor_reduced
    MinorAxisTreduced = MinorAxis_T - minor_reduced
    phi = np.arctan2(y, x)
    rad = np.sqrt(x ** 2 + y ** 2)
    R_prime = np.sqrt((MajorAxisLength * np.cos(phi)) ** 2 
                      + (MinorAxisLength * np.sin(phi)) ** 2)

    scaled_r = rad / R_prime

    scaled_thresh = 0.7* MinorAxisTreduced / MinorAxisLength # For gauss function
    return (scaled_r,scaled_thresh)

def ang_factor(x, y, periods=4*N_gyri-2):
    """
    Computes the radial sine function for the given Cartesian coordinates.

    Parameters:
    - x: X-coordinate (array-like).
    - y: Y-coordinate (array-like).
    - periods: Number of periods of the sine wave.

    Returns:
    - Values of the radial sine function for the given Cartesian coordinates.
    """
    theta = np.arctan2(y, x)  # Compute theta from Cartesian coordinates
    return np.sin(periods * theta)  # Sine function of theta

def growth_rate_total(x, y, thresh, flag=0):
    (scaled_r, scaled_thresh) = scaled_xy(x, y)
    
    # Compute the growth value for r >= threshold
    if flag == 'phi':
        growth_above = 0.5 * gaussian(scaled_r, scaled_thresh, 0.4) * (ang_factor(x, y) + 1.0)
    elif flag == 'H':
        growth_above = 0.5 * heaviside(scaled_r, scaled_thresh, gamma=10) * (ang_factor(x, y) + 1.0)
    else:
        growth_above = 0.0
    
    # Value when r < threshold (constant 1.0)
    growth_below = 1.0
    
    # Apply smooth transition instead of sharp threshold
    gw = apply_smooth_transition(scaled_r, thresh, growth_below, growth_above)
    
    return gw, scaled_r

# Optional: Function to compare sharp vs smooth transitions
def compare_transitions():
    """Generate comparison plot of sharp vs smooth transitions"""
    r_test = np.linspace(0, 2, 1000)
    threshold_test = 0.8
    
    # Sharp transition
    sharp = np.where(r_test < threshold_test, 1.0, 0.0)
    
    # Smooth transitions
    sigmoid = sigmoid_transition(r_test, threshold_test, SMOOTHING_STEEPNESS)
    sigmoid_values = 1.0 * (1 - sigmoid) + 0.0 * sigmoid
    
    tanh = tanh_transition(r_test, threshold_test, SMOOTHING_STEEPNESS)
    tanh_values = 1.0 * (1 - tanh) + 0.0 * tanh
    
    smoothstep = smoothstep_transition(r_test, threshold_test, TRANSITION_WIDTH)
    smoothstep_values = 1.0 * (1 - smoothstep) + 0.0 * smoothstep
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(r_test, sharp, 'k--', label='Sharp threshold', linewidth=3)
    ax.plot(r_test, sigmoid_values, 'b-', label=f'Sigmoid (steepness={SMOOTHING_STEEPNESS})', linewidth=3)
    ax.plot(r_test, tanh_values, 'g-', label=f'Tanh (steepness={SMOOTHING_STEEPNESS})', linewidth=3)
    ax.plot(r_test, smoothstep_values, 'r-', label=f'Smoothstep (width={TRANSITION_WIDTH})', linewidth=3)
    ax.axvline(threshold_test, color='gray', linestyle=':', alpha=0.5, linewidth=2)
    ax.set_xlabel(r'$\tilde{r}$')
    ax.set_ylabel('Transition value')
    ax.set_title(f'Comparison of Smoothing Methods\nSMOOTHING_TYPE = {SMOOTHING_TYPE}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('transition_comparison.pdf', format='pdf', dpi=300)
    plt.show()

# Uncomment to see comparison plot
# compare_transitions()

x_W = np.linspace(0, MajorAxis_W, 500)
y_W = np.linspace(0, MinorAxis_W, 500)
x_G = np.linspace(0, MajorAxis_G, 498)
y_G = np.linspace(0, MinorAxis_G, 498) 
x_G = np.append(x_G,[MajorAxis_W,MinorAxis_W])
x_G.sort()

y_G = np.append(y_G,[MajorAxis_W,MinorAxis_W])
y_G.sort()

(scaled_coord,scaled_thresh) = scaled_xy(x_W, y_W)
theta = np.linspace(0,np.pi/2,num=500) 
periods = 4*N_gyri - 2
angular_f2 = np.sin(periods * theta)

X, Y = np.meshgrid(x_W, y_W)
X_G,Y_G = np.meshgrid(x_G, y_G) 

ellipse = ((X / MajorAxis_W) ** 2) + ((Y / MinorAxis_W) ** 2)
ellipse_G = ((X_G / MajorAxis_G) ** 2) + ((Y_G / MinorAxis_G) ** 2)
ellipse_G_W_val = np.sqrt(((MajorAxis_W / MajorAxis_G) ** 2)) # + ((MinorAxis_W / MinorAxis_G) ** 2))

major_axis_reduced = MajorAxis_W - major_reduced
minor_axis_reduced = MinorAxis_W - minor_reduced
ellipse_reduced = ((X / major_axis_reduced) ** 2) + ((Y / minor_axis_reduced) ** 2)
inside_ellipse = ellipse <= 1
inside_ellipse_gray_strip_W = (ellipse<=1) & (ellipse_G>=ellipse_G_W_val)
inside_ellipse_reduced = ellipse_reduced <=1
inside_ellipse_gray_strip = (ellipse_G<=1) & (ellipse_G>=ellipse_G_W_val)
inside_ellipse_gray = (ellipse_G<=1)

Z = np.full_like(X, np.nan)  # Initialize Z with NaN for outside points
PHI = np.full_like(X, np.nan)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        if inside_ellipse[i, j]:
            Z[i, j], PHI[i, j] = growth_rate_total(X[i, j], Y[i, j], thresh, flag=flag)
        # else:
        #     Z[i,j] = -0.1

viridis = mpl.colormaps['plasma']


# save_dir_figures = "figures"
# os.makedirs(save_dir_figures, exist_ok=True)  # Create directory if it doesn't exist        

xmax_val = np.nanmax(PHI)
yticks = [0,0.5,1.0]

fig, ax = plt.subplots(1,2,figsize=(25,10))

if flag=='phi':
    ### f1 -> USE FOR PHI
    ax[0].plot(scaled_coord,gaussian(scaled_coord, scaled_thresh,0.4),color='black',lw=lw)
    ax[0].set_ylabel(r'$f_1^{\phi}$',labelpad=lp)
else: ### f1 -> USE FOR H
    ax[0].plot(scaled_coord,heaviside(scaled_coord, scaled_thresh,gamma=10),color='black',lw=lw)
    ax[0].set_ylabel(r'$f_1^{{H}}$',labelpad=lp)

ax[0].set_xlim(xmin=0.0,xmax=xmax_val)
ax[0].set_xbound(lower=0.0,upper=xmax_val)
ax[0].set_xlabel(r'$\tilde{r}$ [-]',labelpad=lp)
ax[0].set_yticks(yticks)


##  f2
# Set x-ticks at multiples of pi/4
ax[1].xaxis.set_major_locator(MultipleLocator(base=np.pi / 4))
# Format the x-tick labels as multiples of pi
def format_func(value, tick_number):
    n = int(np.round(value / (np.pi / 4)))
    if n == 0:
        return r"$0$"
    elif n == 1:
        return r"$\frac{\pi}{4}$"
    elif n == 2:
        return r"$\frac{\pi}{2}$"
    elif n == 3:
        return r"$\frac{3\pi}{4}$"
    elif n == 4:
        return r"$\pi$"
    else:
        return r"${}\pi/4$".format(n)

ax[1].xaxis.set_major_formatter(FuncFormatter(format_func))
ax[1].plot(theta,0.5*angular_f2+0.5,color='black',lw=lw)
ax[1].set_yticks(yticks)
ax[1].set_xlabel(r'$\psi$ [rad]',labelpad=lp)
ax[1].set_ylabel(r'$f_2$',labelpad=lp)

plt.tight_layout()
###
fig, ax = plt.subplots(figsize=(8,6))
## f1 * f2
CS = ax.contourf(X, Y, Z, levels=100, cmap='plasma')
ax.contour(X, Y, ellipse, levels=[1], colors='black')

ax.set_ylim(0,MinorAxis_W)
ax.set_xlim(0,MajorAxis_W)
ax.axis('off')
ax.set_aspect('equal')


# Remove ticks and labels
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
plt.tight_layout()
plt.savefig(fig_filename,format='pdf',dpi=dpi)
#%%
PHI_1 = PHI[PHI<=1.0]
# Create the filled contour plot for \tilde{r}
fig, ax = plt.subplots(figsize=(8,6))
# Filled contour
cf = ax.contourf(X, Y, PHI, levels=100, cmap='plasma')

# Outer ellipse (black solid)
ax.contour(X, Y, ellipse, levels=[1], colors='black')

# Add the PHI=1.0 contour as dotted gray
ax.contour(X, Y, PHI, levels=[1.0],
           colors='gray', linestyles='dotted', linewidths=4)

ax.set_ylim(0,MinorAxis_W+1.0)
ax.set_xlim(0,MajorAxis_W+1.0)
ax.axis('off')
ax.set_aspect('equal')
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
plt.savefig('Fig2_A.pdf',format='pdf',dpi=dpi)

plt.tight_layout()
#%%
x_G = np.linspace(0, MajorAxis_G, 498)
y_G = np.linspace(0, MinorAxis_G, 498) 
x_G = np.append(x_G,[MajorAxis_W,MinorAxis_W])
x_G.sort()
y_G = np.append(y_G,[MajorAxis_W,MinorAxis_W])#
y_G.sort()

X_G,Y_G = np.meshgrid(x_G, y_G) 

ellipse_G = ((X_G / MajorAxis_G) ** 2) + ((Y_G / MinorAxis_G) ** 2)
ellipse_G_W_val = np.sqrt(((MajorAxis_W / MajorAxis_G) ** 2)) # + ((MinorAxis_W / MinorAxis_G) ** 2))


inside_ellipse_white = ellipse_G <= ellipse_G_W_val # Points in ellipse that belong to the white matter.
inside_ellipse_gray_strip_W = (ellipse_G<=1) & (ellipse_G>=ellipse_G_W_val) # Points in ellipse that belog to gray matter

inside_ellipse_gray = (ellipse_G<=1) # All points in ellipse

Z_W = np.full_like(X_G, np.nan)  # Initialize Z with NaN for outside points
# PHI = np.full_like(X, np.nan)
for i in range(X_G.shape[0]):
    for j in range(X_G.shape[1]):
        if inside_ellipse_white[i, j]:
            Z_W[i, j], _ = growth_rate_total(X_G[i, j], Y_G[i, j], thresh, flag=flag)
        elif inside_ellipse_gray_strip_W[i,j] and flag=='phi':
            Z_W[i,j] = 0.0
        elif inside_ellipse_gray_strip_W[i,j] and flag=='H':
            Z_W[i,j] = 1.0            
fig, ax = plt.subplots(figsize=(8,6))

# Outer ellipse (black solid)
ax.contourf(X_G, Y_G, Z_W, levels=100,cmap='plasma')
ax.contour(X_G, Y_G, ellipse_G, levels=[ellipse_G_W_val,1.0],colors='gray',linewidths=4)

ax.set_ylim(0,MinorAxis_G+1.0)
ax.set_xlim(0,MajorAxis_G+1.0)
ax.axis('off')
ax.set_aspect('equal')
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# plt.savefig(os.path.join(save_dir_figures, 'GZ_growing.pdf'),format='pdf',dpi=dpi)
plt.savefig('Fig1_phase3_push.pdf',format='pdf',dpi=dpi)

Z_G = np.full_like(X_G, np.nan)  # Initialize Z with NaN for outside points
# PHI = np.full_like(X, np.nan)
for i in range(X_G.shape[0]):
    for j in range(X_G.shape[1]):
        if inside_ellipse_gray_strip[i, j]:
            Z_G[i, j] = 1.0
        elif inside_ellipse_gray[i,j]:
            Z_G[i,j] = 0.0
        # else:
            
fig, ax = plt.subplots(figsize=(8,6))

# Outer ellipse (black solid)
ax.contourf(X_G, Y_G, Z_G, levels=100,cmap='plasma')
ax.contour(X_G, Y_G, ellipse_G, levels=[ellipse_G_W_val,1.0],colors='gray',linewidths=4)

ax.set_ylim(0,MinorAxis_G+1.0)
ax.set_xlim(0,MajorAxis_G+1.0)
ax.axis('off')
ax.set_aspect('equal')
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
# plt.savefig(os.path.join(save_dir_figures, 'cortex_growing.pdf'),format='pdf',dpi=dpi)
plt.savefig('Fig1_phase1_2_cortex_growing.pdf',format='pdf',dpi=dpi)

# Print smoothing parameters used
print(f"\nSmoothing parameters:")
print(f"SMOOTHING_TYPE: {SMOOTHING_TYPE}")
print(f"SMOOTHING_STEEPNESS: {SMOOTHING_STEEPNESS}")
if SMOOTHING_TYPE == 'smoothstep':
    print(f"TRANSITION_WIDTH: {TRANSITION_WIDTH}")
print("\nTo change smoothing behavior, modify the USER-DEFINED SMOOTHING PARAMETERS section at the top of the script.")