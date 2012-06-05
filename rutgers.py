from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.spatial import KDTree, cKDTree
import matplotlib.cm as cm

import pycdf

import base
import projmaps
import gmtgrid
import figpref
from hitta import GBRY

class NWA(base.Njord):
    """Setup North West Atlantic Grid"""
    def __init__(self, datadir="/Volumes/keronHD4/rutgersNWA/",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        super(NWA, self).__init__()
        self.i1 = 0
        self.i2 = 721
        self.j1 = 0
        self.j2 = 361
        self.datadir = datadir
        g = pycdf.CDF(datadir + '/NWA_grd.nc')
        self.llat = g.var('lat_rho')[:]
        self.llon = g.var('lon_rho')[:]-360
        self.depth = g.var('h')[:]
        self.Cs_r = g.var('Cs_r')[:]
        self.region = "nwa_small"
        self.add_mp()

    def load(self,fldname,jd=730217,yr=0,mn=1,dy=1,hr=3):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy,hr))
        filename = "/nwa_avg_%05i.nc" % (jd - 714782)
        nc = pycdf.CDF("%s/%s/%s" % (self.datadir, pl.num2date(jd).year,
                                     filename))     
        fld =  np.squeeze(nc.var(fldname)[:])
        fld[fld>9999] = np.nan
        self.__dict__[fldname] = fld
        self.ssh =  np.squeeze(nc.var('zeta')[:])
        self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
                     self.Cs_r[:,np.newaxis,np.newaxis])

    def add_landmask(self):
        g = pycdf.CDF(self.datadir + '/NWA_grd.nc')
        self.landmask = g.var('mask_rho')[:]

    def add_utmxy(self):
        g = pycdf.CDF(self.datadir + '/NWA_grd.nc')
        self.utmx = g.var('x_rho')[:]
        self.utmy = g.var('y_rho')[:]

class Coral(base.Njord):
    """Setup INdonesial flowthrough Grid"""
    def __init__(self, datadir="/projData/rutgers/CORAL/",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        super(Coral, self).__init__()
        self.i1 = 0
        self.i2 = 1281
        self.j1 = 0
        self.j2 = 641
        self.datadir = datadir
        g = pycdf.CDF(datadir + '/coral_grd.nc')
        self.llat = g.var('lat_rho')[:]
        self.llon = g.var('lon_rho')[:]
        self.depth = g.var('h')[:]
        self.Cs_r = g.var('Cs_r')[:]
        self.region = "indthr"
        self.add_mp()

    def load(self,fldname,jd=730217,yr=0,mn=1,dy=1,hr=3):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy,hr))
        filename = "/nwa_avg_%05i.nc" % (jd - 714782)
        nc = pycdf.CDF("%s/%s/%s" % (self.datadir, pl.num2date(jd).year,
                                     filename))     
        fld =  np.squeeze(nc.var(fldname)[:])
        fld[fld>9999] = np.nan
        self.__dict__[fldname] = fld
        self.ssh =  np.squeeze(nc.var('zeta')[:])
        self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
                     self.Cs_r[:,np.newaxis,np.newaxis])

    def add_landmask(self):
        g = pycdf.CDF(self.datadir + '/coral_grd.nc')
        self.landmask = g.var('mask_rho')[:]

    def add_utmxy(self):
        g = pycdf.CDF(self.datadir + '/coral_grd.nc')
        self.utmx = g.var('x_rho')[:]
        self.utmy = g.var('y_rho')[:]
        

