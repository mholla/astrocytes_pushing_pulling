#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 15:33:11 2025

@author: ktaneja
"""
import numpy as np
from scipy.spatial import cKDTree

def average_shortest_distance(array1, array2):
    """Average shortest distance from array1 to array2 using KD-Tree (fast).
    This is an optimized way for finding the closest distance from B to A.
    """
    A = np.asarray(array1, float)
    B = np.asarray(array2, float)
    d, idx = cKDTree(B).query(A, k=1, workers=-1)
    neighbors = B[idx]
    return float(np.mean(d)), d, idx, neighbors


def polyline_length(points):
    """Sum of consecutive segment lengths for an open polyline."""
    pts = np.asarray(points, float)
    if len(pts) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())


def point_to_polyline_min_dists(points, polyline, eps=1e-15):
    """
    Vectorized minimum distance from each point to an open polyline.
    Given a bunch of points, goal is to find the minimum euclidean distance (straight line)
    to the polyline (a bunch of connected straight line segments). 
    A point P_i is projected on every segment A_j,B_j of the polyline, 
    then it's relative position on A_j,B_j is found and then the distance 
    between the relative position and P_i is calculated.   
    

    Returns
    -------
    (N,) array of min distances for each point in `points`.
    """
    P = np.asarray(points, float)
    C = np.asarray(polyline, float)


    A, B = C[:-1], C[1:]              # (M-1,2) each: segment endpoints
    AB = B - A                        # (M-1,2): segment vectors
    AB2 = np.einsum('ij,ij->i', AB, AB)  # (M-1,): squared lengths ||AB||^2 = AB2_i = AB_ik AB_ik
    
    AP = P[:, None, :] - A[None, :, :] # (N, M-1, 2): P_i - A_j via broadcasting
    denom = np.where(AB2 < eps, 1.0, AB2) #np.where(condition Where True yield x otherwise yield y, x, y)
    t = np.einsum('nij,ij->ni', AP, AB) / denom    # (N, M-1) -> Length of projection of vector AP on AB, scaled by length of vector AB, i.e., relative position.
    t = np.where(AB2 < eps, 0.0, t)                # if degenerate, force t=0 (use A_j)
    t = np.clip(t, 0.0, 1.0)                       # clamp to the segment

    closest = A[None, :, :] + t[:, :, None] * AB[None, :, :] # (N, M-1, 2), position of closest point to P on AB
    dists = np.linalg.norm(P[:, None, :] - closest, axis=2) # Distance
    return dists.min(axis=1)


def point_to_polyline_min_dists_with_foot(points, polyline, eps=1e-15):
    """
    Vectorized min distance from each point to an open polyline + foot of perpendicular.
    Given a bunch of points, goal is to find the minimum euclidean distance (straight line)
    to the polyline (a bunch of connected straight line segments). 
    A point P_i is projected on every segment A_j,B_j of the polyline, 
    then it's relative position on A_j,B_j is found and then the distance 
    between the relative position and P_i is calculated.   
    
    Returns
    -------
    dmin : (N,) array
    foot : (N,2) array
    seg_idx : (N,) int array
    """
    P = np.asarray(points, float)
    C = np.asarray(polyline, float)

    A, B = C[:-1], C[1:]              # (M-1,2) each: segment endpoints
    AB = B - A                        # (M-1,2): segment vectors
    AB2 = np.einsum('ij,ij->i', AB, AB)  # (M-1,): squared lengths ||AB||^2 = AB2_i = AB_ik AB_ik
    
    AP = P[:, None, :] - A[None, :, :] # (N, M-1, 2): P_i - A_j via broadcasting
    denom = np.where(AB2 < eps, 1.0, AB2) #np.where(condition Where True yield x otherwise yield y, x, y)
    t = np.einsum('nij,ij->ni', AP, AB) / denom    # (N, M-1) -> Length of projection of vector AP on AB, scaled by length of vector AB, i.e., relative position.
    t = np.where(AB2 < eps, 0.0, t)                # if degenerate, force t=0 (use A_j)
    t = np.clip(t, 0.0, 1.0)                       # clamp to the segment

    closest = A[None, :, :] + t[:, :, None] * AB[None, :, :] # (N, M-1, 2), position of closest point to P on AB
    dists = np.linalg.norm(P[:, None, :] - closest, axis=2) # Straight line distance


    # Index, for each point P_i (row), of the segment j that gives the MIN distance.
    # dists has shape (N, M-1), so argmin over axis=1 returns an array (N,)
    # where jmin[i] ∈ {0, ..., M-2} is the index of the best segment for P_i.
    jmin = np.argmin(dists, axis=1)

    # Row indices [0, 1, ..., N-1] used to gather the per-point minima cleanly.
    rows = np.arange(P.shape[0])

    # Extract the minimum distance for each point using (row i, best-segment jmin[i]).
    # Result shape: (N,)
    dmin = dists[rows, jmin]

    # Extract the corresponding closest-point (foot of perpendicular) on the polyline
    # for each P_i, again indexing with (row i, best-segment jmin[i], :)
    # Result shape: (N, 2)
    foot = closest[rows, jmin, :]

    # Keep the best segment index as int; seg_idx[i] tells which segment [A_j, B_j]
    # contains the foot for point P_i. Useful if you need the segment later.
    seg_idx = jmin.astype(int)

    # Return:
    # - dmin: per-point minimum distances to the polyline        (N,)
    # - foot: per-point closest coordinates on the polyline      (N, 2)
    # - seg_idx: per-point segment indices where the foot lies   (N,)
    return dmin, foot, seg_idx
