
import numpy as np
import pylab as pl
import pandas as pl

import lldist

def approx_dx(lonarr, latarr):
    """Caclulate dx from lon and lat"""
    if lonarr.ndim == latarr.ndim == 1:
        lonarr = lonarr[None, :]
        latarr = lonarr[:, None]
        
    dx = lonarr * np.nan
    for i in range(lonarr.shape[0]):
        latvec1 = (latarr[i,0:-2] + latarr[i,1:-1]) / 2
        latvec2 = (latarr[i,2:]   + latarr[i,1:-1]) / 2
        lonvec1 = (lonarr[i,0:-2] + lonarr[i,1:-1]) / 2
        lonvec2 = (lonarr[i,2:]   + lonarr[i,1:-1]) / 2
        dx[i,1:-1] = lldist.ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2)
    dx[:, 0] = 2 * dx[:,1]  - dx[:,2]
    dx[:,-1] = 2 * dx[:,-2] - dx[:,-3]
    return dx

def approx_dy(lonarr, latarr):
    """Caclulate dy from lon and lat"""
    if lonarr.ndim == latarr.ndim == 1:
        lonarr = lonarr[None, :]
        latarr = lonarr[:, None]
    
    dy = latarr * np.nan
    for j in range(lonarr.shape[1]):
        latvec1 = (latarr[0:-2,j] + latarr[1:-1,j]) / 2
        latvec2 = (latarr[2:,j]   + latarr[1:-1,j]) / 2
        lonvec1 = (lonarr[0:-2,j] + lonarr[1:-1,j]) / 2
        lonvec2 = (lonarr[2:,j]   + lonarr[1:-1,j]) / 2
        dy[1:-1,j] = lldist.ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2)
    dy[ 0,:] = 2 * dy[ 1,:] - dy[ 2,:]
    dy[-1,:] = 2 * dy[-2,:] - dy[-3,:]
    return dy
