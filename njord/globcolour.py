import os, os.path
import glob
import subprocess as sbp
from datetime import datetime as dtm
import urllib
import urllib2
import re
from distutils import spawn
import warnings

import numpy as np
import pylab as pl

import netCDF4

import base
from njord.utils import yrday

class Base(base.Grid):

    def __init__(self, **kwargs):
        super(Base, self).__init__(**kwargs)

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        with netCDF4.Dataset(self.gridfile) as nc:
            self.latvec = nc.variables["lat"][:]
            self.lonvec = nc.variables["lon"][:]
            self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    @property
    def filedict(self):
        if not hasattr(self, "_filedict"):
            filestamp = "L3m_*-*__GLOB_4_%s-*_CHL1_8D_00.nc" % self.algo
            flist = glob.glob(os.path.join(self.datadir, filestamp))
            self._filedict  = {}
            for fn in flist:
                filename = os.path.basename(fn)
                self._filedict[pl.datestr2num(filename[4:12])] = fn
        return self._filedict

    @property
    def fulltvec(self):
        return np.sort(np.array(self.filedict.keys()))

    @property
    def minjd(self):
        return self.fulltvec.min()
    
    @property
    def maxjd(self):
        return self.fulltvec.max()
            
    def load(self, fld="chl", fldtype="8D", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f and field=%s" % (self.jd, fld))
        self.filename = self.filedict[
            self.fulltvec[self.fulltvec <= self.jd].max()]
        with netCDF4.Dataset(os.path.join(self.datadir, self.filename)) as nc:
            nc.set_auto_mask(False)
            fld = nc.variables["CHL1_mean"][self.j1:self.j2,self.i1:self.i2].copy()
            fld[fld<0] = np.nan
        setattr(self, "chl", fld)
        self.vprint( "datadir:        %s" % self.datadir)
        self.vprint( "filename:       %s" % os.path.basename(self.filename))

    @property
    def landmask(self):
        """Return a landmask"""
        if not hasattr(self,'_landmask'):
            try:
                fld = self.get_field('par', fldtype='CU', nan="nan")
            except KeyError:
                fld = self.get_field('chl', fldtype='CU', nan="nan")
            self._landmask = np.isnan(fld)
        return self._landmask

    @property
    def fieldlist(self):
        return ["chl",]

    
class AVW(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.algo = "AVW"
        super(AVW, self).__init__(**kwargs)
