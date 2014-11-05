import os, os.path
import urllib
import zipfile
import tempfile

import numpy as np

import base
from utils import yrday

class Glob30s(base.Grid):

    def __init__(self, **kwargs):
        super(Glob30s, self).__init__(**kwargs)
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
    
    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        dllarr = lambda dll,mt: np.cumsum([dll] * mt)
        self.lonvec = (self.full_lon1 + dllarr(self.dlon, self.imt))
        self.latvec = (self.full_lat2 - dllarr(self.dlat, self.jmt))
        self._llboundaries_to_ij()
        self.full_imt = self.imt
        self.full_jmt = self.jmt
        self.imt = self.i2 - self.i1 + 1
        self.jmt = self.j2 - self.j1 + 1
                
    def load(self, fld, mn=None, **kwargs):
        """Load the data field associated with a given time."""
        self._timeparams(**kwargs)
        self.mn = mn if mn is not None else self.mn
        if fld == "alt":
            fn = "%s/alt.bil" % self.datadir
        else:
            fn = "%s/tmean_%i.bil" % (self.datadir, self.mn)
        if not os.path.isfile(fn):
            self.download(fld)
        mat = np.fromfile(fn, dtype="<i2").reshape(self.full_jmt, self.full_imt)
        mat = mat[self.j1:self.j2, self.i1:self.i2].astype(np.float16)
        mat[mat < -9998] = np.nan
        setattr(self, fld, mat * self.vc['alt'][2])
        
    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'landmask'):
            self.load('par','CU', nan="nan")
            self.landmask = np.isnan(self.par)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return  {'tmean': ['Mean Temperature', 'C', 0.1],
                 'tmax':  ['Max Temperature',  'C', 0.1],
                 'prec':  ['Precipitation',   'mm', 1],
                 'alt':   ['Altitide',         'm', 1],
                   }

    def download(self, fldname):
        """Download missing data files"""
        print "Downloading all %s datafiles. This will take a while." % fldname
        url = "%s/%s_30s_bil.zip" % (self.dataurl, fldname)
        tmpzipfile = os.path.join(tempfile.gettempdir(), "tmp.zip")
        urllib.urlretrieve(url, tmpzipfile)
        with zipfile.ZipFile(tmpzipfile, "r") as z:
            z.extractall(self.datadir)

    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            self._landmask =  ~np.isnan(self.get_field('alt'))
        return self._landmask        

  

def cristian(filename="growth_days.csv", tbase=0):
    import pandas as pd
    from scipy.spatial import cKDTree

    
    ll = pd.read_csv('Coordinates_SLU.csv')
    wc = Glob30s(lat1=53,lat2=73,lon1=0,lon2=30)
    kd = cKDTree(np.vstack((wc.llon.flat,wc.llat.flat)).T)
    dist,ij = kd.query(np.vstack((ll['Xswr99'], ll['Yswr99'])).T)
    ijmask = ij < wc.llat.size
    mnln = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    temp = np.zeros((12,len(ij[ij<wc.llon.size])))
    for n,mn in enumerate(np.arange(1,13)):
        print mn
        wc.load('tmean', mn=mn)
        temp[n,:] = wc.tmean.flat[ij[ijmask]]
    grday = temp.copy()
    grday[grday<tbase] = np.nan
    grdays = np.nansum((grday-tbase) *  mnln[:,np.newaxis],axis=0)
    
    mat = np.vstack(( ll['Xswr99'][ijmask], ll['Yswr99'][ijmask], grdays))
    headerstr = ("Lon,Lat,Growing_days")
    np.savetxt(filename, mat.T, fmt="%.8f", delimiter=",", header=headerstr) 

    return ll['Xswr99'], ll['Yswr99'], grdays


"""
from scipy.interpolate import Rbf, SmoothBivariateSpline
from njord import worldclim

wc = worldclim.Glob30s(lat1=53,lat2=73,lon1=3,lon2=43)
mat = loadtxt('S_lakes_Bror.csv', delimiter=",",skiprows=1,usecols=[1,2,3,4,5,6])
 rbf = Rbf(mat[1:500,0],mat[1:500,1],mat[1:500,4], epsilon=2)
 ZI = rbf(wc.llon, wc.llat))
"""
