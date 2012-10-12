from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
import matplotlib.cm as cm
from scipy.io import netcdf_file, netcdf_variable

import base

class Glob_025_ll(base.Grid):
    """Setup Ecco 0.25 deg fields."""
    def __init__(self, **kwargs):
        super(Glob_025_ll, self).__init__()
        self.add_mp()

    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile, 'r')
        self.lat = g.variables['LATITUDE_T'][:]
        self.lon = g.variables['LONGITUDE_T'][:]-self.lonoffs
        self.depth = g.variables['DEPTH_T'][:]
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.vclist = {'u':'UVEL','v':'VVEL','w':'WVEL'}

    def load(self,fldname,jd=732170,yr=0,mn=1,dy=1,hr=3):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy,hr))
        jd = np.round((732170.0+1)/3)*3-1
        dtm = pl.num2date(jd)
        filename = ("%s.1440x720x50.%04i%02i%02i.nc" %
                    (self.vclist[fldname],dtm.year,dtm.month,dtm.day))
        nc = netcdf_file("%s/%s" % (self.datadir, filename))     

        fld =  np.squeeze(nc.variables[self.vclist[fldname]][:]).copy()
        fld[fld==0] = np.nan
        self.__dict__[fldname] = fld
        #self.ssh =  np.squeeze(nc.var('zeta')[:])
        #self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
        #             self.Cs_r[:,np.newaxis,np.newaxis])

    def add_landmask(self):
        """ Generate a landmask for t-positions"""
        self.landmask = self.depth.copy()
    
