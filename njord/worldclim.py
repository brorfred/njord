import os, os.path
from datetime import datetime as dtm
import urllib
import urllib2
import re

import numpy as np
import pylab as pl

import base
import yrday

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
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.mn = mn if mn is not None else mn
        print self.mn
        fn = "%s/tmean_30s_bil/tmean_%i.bil" % (self.datadir, self.mn)
        mat = np.fromfile(fn, dtype="<i2").reshape(self.full_jmt, self.full_imt)
        mat = mat[self.j1:self.j2, self.i1:self.i2].astype(np.float16)
        mat[mat < -9998] = np.nan
        setattr(self, fld, mat/10)
        
    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'landmask'):
            self.load('par','CU', nan="nan")
            self.landmask = np.isnan(self.par)
        
    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['KD490_Kd_490',         'km', 0,      100],
                   'chl': ['CHL_chlor_a',          'km', 0.001,  200],
                   'poc': ['POC_poc',              'km', 0,    10000],
                   'bbp': ['GIOP01_bbp_443_giop',  'km', 0,    10000],
                   'sst': ['SST',                  '',  -2,       45],
                   'par': ['PAR_par',              'km',-5,      200],
                   'flh': ['FLH_nflh',             'km', 0,    10000],
                   'ipar': ['FLH_ipar',            'km', 0,    10000],
                   'eup': ['KDLEE_Zeu_lee',        'km', 0,    10000],
                   'eup_lee': ['KDLEE_Zeu_lee',    'km', 0,    10000],
                   'eup_mor': ['KDMOREL_Zeu_morel','km', 0,    10000]
                   }

    def download(self, filename):
        """Download a missing file from GSFC's website"""
        uri = "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/"
        url = "%s%s.bz2" % (uri, os.path.basename(filename))
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(filename + ".bz2",'wb')
        output.write(response.read())
        output.close()

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
