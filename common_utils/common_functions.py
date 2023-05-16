# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:42:06 2023

@author: sunney
"""

import numpy as np

def remove_points(ndp_preProcessed):
    
    """ removes dominated points and any duplicates """
    ndp_preProcessed = np.round(ndp_preProcessed, 2)
    NDP = np.unique(ndp_preProcessed, axis=0)
    if ndp_preProcessed[0][0] == -1 and ndp_preProcessed[0][1] == -1:
        NDP = np.delete(NDP, 0, axis=0)
    NDP = keep_efficient(-NDP)
    NDP = -NDP
    return NDP


def keep_efficient(pts):
    'returns Pareto efficient row subset of pts'
    # sort points by decreasing sum of coordinates
    pts = pts[pts.sum(1).argsort()[::-1]]
    # initialize a boolean mask for undominated points
    # to avoid creating copies each iteration
    undominated = np.ones(pts.shape[0], dtype=bool)
    for i in range(pts.shape[0]):
        # process each point in turn
        n = pts.shape[0]
        if i >= n:
            break
        # find all points not dominated by i
        # since points are sorted by coordinate sum
        # i cannot dominate any points in 1,...,i-1
        undominated[i+1:n] = (pts[i+1:] > pts[i]).any(1)  ##
        # keep points undominated so far
        pts = pts[undominated[:n]]
    return pts
