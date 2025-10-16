#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: ktaneja

Quick parameter guide
---------------------
Global plotting / typography
- SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE, LARGE_SIZE: Matplotlib font sizes.
- plt.rc / plt.rcParams: Global fonts (Times New Roman) and TeX rendering. Requires LaTeX in PATH.

IO and behavior flags
- file_path: Tab-delimited XY data (ImageJ export). Selecting this path controls the dataset-specific
             constants set below (thresholds, limits, anchors, etc.).
- file_path_pic: Histology image for overlays (used when pic_flag == 1).
- pic_flag: 1 enables the image overlay figure; 0 skips it.
- fig_filename: Output figure name for the image overlay figure when pic_flag == 1.
- save_dir_figures: Output directory for figures (created if missing).
- dpi: Figure save resolution (used at the very end when saving overlay figure).

Dataset-specific tuning (auto-selected by file_path)
- ylim: Upper Y cutoff for filtering input data from image.
- y_thresh, y_thresh_left, y_thresh_right: Thresholds used to remove rim points on left/right ends.
- xmin_pic, ymin_pic: Lower-left anchor for placing the histology image in plot coordinates.
- y_thresh_pic, x_thresh_pic: Extent offsets to align the histology image with data coordinates.
- alphas: List of smoothing parameters (s) for scipy.interpolate.splprep; larger s -> smoother curve.
- pt1_x, pt1_y / pt2_x, pt2_y: Crown reference points used to find start/end indices on the spline.
- s_size: Scatter marker size for the raw CP points in spline previews.
- (Optional) dum_ptx, dum_pty: Extra anchors (left in place for potential future use).

Metrics & outputs produced
- GI: Gyrification index proxy = length(pial curve) / length(convex hull of pial curve).

"""


# -------------------------------
# 1) Imports & global style
# -------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree
from scipy.interpolate import splprep, splev
from matplotlib.collections import LineCollection
import os
from pathlib import Path
from Gyrification_functions import *
# Make sure LaTeX binaries are on PATH if using usetex=True
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'  # '/latex'

# Global typography sizes
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 24
LARGE_SIZE = 18

# Matplotlib rcParams for consistent sizing/typography across all figures
plt.rc('font', size=BIGGER_SIZE)           # default text
plt.rc('axes', titlesize=BIGGER_SIZE)      # axes title
plt.rc('axes', labelsize=BIGGER_SIZE)      # axes labels
plt.rc('xtick', labelsize=BIGGER_SIZE)     # x tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)     # y tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)     # legend text
plt.rc('figure', titlesize=LARGE_SIZE)     # figure title
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('text', usetex=True)                # requires system LaTeX install

# -------------------------------
# 2) Raw Image Loading
# -------------------------------
current_folder = Path(__file__).resolve().parent
main_directory = current_folder.parent
image_data_folder = main_directory/"Raw_extracted_points"
raw_image_folder = main_directory/"Raw_images"

file_path = image_data_folder/"P16_2_XYdata_CP_raw.txt"  # Input CP coordinates (tab-delimited from ImageJ)
pic_flag = 1                                   # 1=overlay hull/surfaces on histology image
file_path_pic = raw_image_folder/'P16_2.png'   # Histological image name

fig_filename = 'P16_2_hull.pdf'  # Output filename when pic_flag==1


save_dir_figures = "figures"                       # Figure output directory
os.makedirs(save_dir_figures, exist_ok=True)
#

dpi = 600*2*2  # High-resolution output 


# -------------------------------
# 3) Dataset-specific tuning 
# -------------------------------
if file_path == image_data_folder/"P16_1_XYdata_CP_raw.txt":
    ylim = 12.6
    y_thresh = 0.75
    y_thresh_left = y_thresh
    y_thresh_right = y_thresh
    xmin_pic = 0.0
    ymin_pic = 0.0
    y_thresh_pic = 0.45
    x_thresh_pic = -0.5
    alphas = [100, 50]    # spline smoothing levels to preview, uses 50 eventually
    pt1_x = 1.3854
    pt1_y = ylim
    pt2_x = 19.9688
    pt2_y = ylim
    s_size = 0.01

if file_path == image_data_folder/"P16_2_XYdata_CP_raw.txt":
    ylim = 13.4 + 0.5 + 0.5
    y_thresh = 0.5
    y_thresh_left = y_thresh + 0.5
    y_thresh_right = y_thresh
    xmin_pic = 0.4
    ymin_pic = 0.0
    y_thresh_pic = -0.45
    x_thresh_pic = 0.25
    alphas = [100, 50]
    pt1_x = 0.8
    pt1_y = ylim
    pt2_x = 14.5 + 5.0
    pt2_y = ylim
    s_size = 0.01

if file_path == image_data_folder/"P16_3_XYdata_CP_raw.txt":
    ylim = 13.4
    y_thresh = 0.25
    y_thresh_left = y_thresh + 1.0
    y_thresh_right = y_thresh + 0.5
    xmin_pic = 0.0
    ymin_pic = 0.0
    # y_thresh_pic = -0.95
    # x_thresh_pic = 0.1#
    y_thresh_pic = -0.2
    x_thresh_pic = -0.2#
    alphas = [100, 50]
    pt1_x = 3.0 - 1.5
    pt1_y = ylim
    pt2_x = 17.5 + 1.5
    pt2_y = ylim
    s_size = 0.01

if file_path == image_data_folder/"P10_1_XYdata_CP_raw.txt":
    ylim = 10.6
    y_thresh = 0.25
    y_thresh_left = y_thresh
    y_thresh_right = y_thresh
    xmin_pic = 0.0
    ymin_pic = 0.0
    y_thresh_pic = y_thresh - 0.45 
    x_thresh_pic = 0.4
    alphas = [25, 30]
    pt1_x = 1.5
    pt1_y = ylim
    pt2_x = 18.5 - 0.5
    pt2_y = ylim
    s_size = 0.01

if file_path == image_data_folder/"P10_2_XYdata_CP_raw.txt":
    ylim = 10.6
    y_thresh = 0.25
    y_thresh_left = y_thresh
    y_thresh_right = y_thresh
    xmin_pic = 0.0
    ymin_pic = 0.0
    y_thresh_pic = y_thresh - 0.5 
    x_thresh_pic = -0.45
    alphas = [25, 30]
    pt1_x = 0.1
    pt1_y = ylim - 1.0
    pt2_x = 20
    pt2_y = ylim
    s_size = 0.01

# -------------------------------
# Geometry helper function to get nearest point within distance
# -------------------------------
def nearest_point_within_radius(center, radius, points):
    """
    Find the nearest point in `points` to `center` within a given `radius`.

    Parameters
    ----------
    center : (2,) array-like
        Target (x,y).
    radius : float
        Radial search bound in data units.
    points : (N,2) array-like
        Candidate points.

    Returns
    -------
    (2,) ndarray or None
        Nearest point to `center` within `radius`, or None if none found.
    """
    points_array = np.array(points)
    distances = np.sqrt(np.sum((points_array - center) ** 2, axis=1))
    within_radius_indices = np.where(distances <= radius)[0]
    if not within_radius_indices.size:
        return None
    nearest_index = within_radius_indices[np.argmin(distances[within_radius_indices])]
    return points_array[nearest_index]

# -------------------------------
# Load & filter CP points; spline preview of cortical plate
# -------------------------------
data = pd.read_csv(file_path, delimiter="\t", header=None, names=["X", "Y"], encoding_errors="ignore")
data = data[data['Y'] < ylim]                  # keep rows below `ylim` cutoff
points = data[['X', 'Y']].values               # (N,2)

# Preview spline for each smoothing level in `alphas`
fig, ax = plt.subplots(len(alphas), figsize=(7, 4 * len(alphas)))
for idx, alpha in enumerate(alphas):
    # Create 2D parametric spline of CP points with smoothing s=alpha
    spl, u = splprep([points[:, 0], points[:, 1]], s=alpha)
    ax[idx].scatter(points[:, 0], points[:, 1], s=s_size, color='blue')
    if spl:
        x, y = splev(u, spl)
        ax[idx].plot(x, y, color='red')

    ax[idx].set_title(f'Alpha = {alpha}')
    ax[idx].set_xlabel('X')
    ax[idx].set_ylabel('Y')

    # Identify left/right crown points on the spline (nearest to pt1/pt2)
    center1 = np.array([pt1_x, pt1_y])
    radius = 1
    points_hull = np.vstack((x.tolist(), y.tolist())).T
    nearest_start = nearest_point_within_radius(center1, radius, points_hull)
    start_ID = np.where((points_hull[:, 0] == nearest_start[0]) & (points_hull[:, 1] == nearest_start[1]))
    ax[idx].plot(nearest_start[0], nearest_start[1], 'y-*')

    center2 = np.array([pt2_x, pt2_y])
    nearest_end = nearest_point_within_radius(center2, radius, points_hull)
    end_ID = np.where((points_hull[:, 0] == nearest_end[0]) & (points_hull[:, 1] == nearest_end[1]))
    ax[idx].plot(nearest_end[0], nearest_end[1], 'g-*')

plt.gca().invert_yaxis()
plt.tight_layout()
#%%
# -------------------------------
# Split spline into GM (pial) segment
# -------------------------------
if end_ID[0][0] < start_ID[0][0]:  # P16_2-like ordering:  GM wraps around

    GM_IDs = np.arange(start_ID[0][0], len(points_hull)).tolist()
    GM_IDs.extend(np.arange(0, end_ID[0][0]).tolist())
    GM_pts = points_hull[GM_IDs]
else:  # P16_1, P16_3-like ordering
    GM_pts = points_hull[start_ID[0][0]:end_ID[0][0], :]

if file_path == image_data_folder/"P10_1_XYdata_CP_raw.txt":
    GM_pts = points_hull[end_ID[0][0]:start_ID[0][0],:]

    
if file_path == image_data_folder/"P10_2_XYdata_CP_raw.txt":
    GM_pts = points_hull[end_ID[0][0]:start_ID[0][0],:]

# Remove rim points by side-specific Y cutoffs (limit crown/edge influence)
GM_mask_left = np.logical_and.reduce([
    GM_pts[:, 1] > (ylim - y_thresh_left),
    GM_pts[:, 0] >= 0,
    GM_pts[:, 0] <= 5,
])
GM_pts = GM_pts[~GM_mask_left]

GM_mask_right = np.logical_and.reduce([
    GM_pts[:, 1] > (ylim - y_thresh_right),
    GM_pts[:, 0] >= 15,
    GM_pts[:, 0] <= 20,
])
GM_pts = GM_pts[~GM_mask_right]


# Quick check plot GM segmentation
plt.figure()
plt.plot(GM_pts[:, 0], GM_pts[:, 1], '-*', label='GM')
plt.gca().invert_yaxis()
plt.legend()



# -------------------------------
# 4) Convex hull of GM (sulcal envelope)
# -------------------------------
gray_surface_hull = ConvexHull(GM_pts)
gray_surface_hull_vertices = gray_surface_hull.vertices
gray_surface_hull_vertices = np.sort(gray_surface_hull_vertices)
gray_surface_hull_coords = GM_pts[gray_surface_hull_vertices]

# Plot style for subsequent figures
ms = 3
lw = 3

# -------------------------------
# 5) Gyrification index (length ratio)
# -------------------------------
L_oc = polyline_length(GM_pts)                  # outer (pial) curve length
L_ch = polyline_length(gray_surface_hull_coords)  # convex hull length
GI = L_oc / L_ch
print('GI=', GI)

#%%
from mpl_toolkits.axes_grid1 import make_axes_locatable
# -------------------------------
# 6) Visualize GM and convex hull
# -------------------------------
plt.figure()
fig, axs = plt.subplots(1, 1, figsize=(10, 10))

axs.plot(GM_pts[:-1, 0], GM_pts[:-1, 1], label='GM', lw=lw, color='lightgrey')


axs.plot(gray_surface_hull_coords[:, 0],
         gray_surface_hull_coords[:, 1],
         '--', color='red', markersize=ms, lw=lw, label='Convex hull')



# cosmetics
axs.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=True)
x_min1, x_max1 = axs.get_xlim()
y_min1, y_max1 = axs.get_ylim()
axs.set_xticks([]); axs.set_yticks([])
axs.axis('off')
axs.set_aspect('equal', adjustable='box')
axs.invert_yaxis()

plt.tight_layout()
plt.show()
#%%

# -------------------------------
# 7) Optional overlay on histology image (pic_flag)
# -------------------------------
if pic_flag == 1:

    img = mpimg.imread(file_path_pic)
    ms = 5
    lw = 3

    fig, axs = plt.subplots(1, 1, figsize=(8, 8))
    # Place image with coordinate extent aligned to data + small x/y offsets
    axs.imshow(img, origin='upper', extent=[xmin_pic, x_max1 + x_thresh_pic, y_max1 + y_thresh_pic, ymin_pic])

    axs.plot(GM_pts[:-1, 0], GM_pts[:-1, 1], color='lightgrey', lw=lw, label='Pial surface')
    axs.plot(gray_surface_hull_coords[:, 0],
             gray_surface_hull_coords[:, 1],
             '--', color='red', markersize=ms, lw=lw, label='Convex hull')


    axs.set_yticks([])
    axs.set_xticks([])
    axs.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1, frameon=True)
    plt.tight_layout()

    # Save high-res overlay figure (PDF)
    plt.savefig(os.path.join(save_dir_figures, fig_filename), format='pdf', dpi=dpi)