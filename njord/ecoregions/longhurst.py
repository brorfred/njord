
import os
import glob

import numpy as  np
import pylab as  pl
import pandas as pd
from scipy.spatial import cKDTree

from njord import rutgers
import ogr2array

def ll2regid(lonvec, latvec, pixel_size=1/12.):
    """Identify Longhurst regions using lat-lon vecs"""
    fn = "data/longhurst_v4_2010/Longhurst_world_v4_2010.shp"
    sh = ogr2array.Shape2Array(fn, pixel_size=pixel_size)
    mat = sh.rasterize().astype(int)
    kd = cKDTree(np.vstack((sh.llon.flat, sh.llat.flat)).T)
    dist,ij = kd.query(np.vstack((lonvec, latvec)).T)
    return mat.flat[ij]

def provcodes():
    """get providence codes"""
    fn = "data/longhurst_v4_2010/Longhurst_world_v4_2010.shp"
    sh = ogr2array.Shape2Array(fn, pixel_size=1/4.)
    return sh.get_fields()["ProvCode"]

def provnames():
    """Get list of long names of all providences"""
    fn = "data/longhurst_v4_2010/Longhurst_world_v4_2010.shp"
    sh = ogr2array.Shape2Array(fn, pixel_size=1/4.)
    return sh.get_fields()["ProvDescr"]

def provsections():
    """Get list of large-scale sections of providences"""
    return [nm.split("-")[0].rstrip() for nm in provnames()]
