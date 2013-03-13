import os, urllib2
from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.io import netcdf_file, netcdf_variable

import requests
import base

class Rutgers(base.Grid):
    """ Baseclass for all rutgers projects """
    def __init__(self, **kwargs):
        super(Rutgers, self).__init__(**kwargs)
        
    def add_landmask(self):
        """ Add a landmask attribute """
        g = netcdf_file(self.gridfile)
        self.landmask = g.variables['mask_rho'][:]

    def add_utmxy(self):
        """Add catersian coordinates """
        g = netcdf_file(self.gridfile)
        self.utmx = g.variables['x_rho'][:]
        self.utmy = g.variables['y_rho'][:]

    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile, 'r')
        self.llat = g.variables['lat_rho'][:]
        self.llon = g.variables['lon_rho'][:]-360
        self.depth = g.variables['h'][:]
        self.Cs_r = g.variables['Cs_r'][:]

    def load(self,fldname,jd=None,yr=0,mn=1,dy=1,hr=3):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy,hr))
        nc = netcdf_file(self.jd2filename(jd))
        fld =  np.squeeze(nc.variables[fldname][:]).copy()
        fld[fld>9999] = np.nan
        self.__dict__[fldname] = fld
        self.ssh =  np.squeeze(nc.variables['zeta'][:])
        self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
                     self.Cs_r[:,np.newaxis,np.newaxis])

class NWA(Rutgers):
    """Setup North West Atlantic"""
    def __init__(self, **kwargs):
        super(NWA, self).__init__(**kwargs)

    def jd2filename(self,jd):
        if jd == None: jd = 730217 
        filename = "/nwa_avg_%05i.nc" % (jd - 714782)
        return "%s/%s/%s" % (self.datadir, pl.num2date(jd).year, filename)    

    def download(self, filename):
        """Download a missing file from Rutger's website"""
        filename = os.path.basename(filename)
        uri = ("http://oceanus.esm.rutgers.edu:8080/" +
               "thredds/fileServer/ROMS/NwA/Run01/Output/Daily/")
        url = "%s%s" % (uri, filename)
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(os.path.join(self.datadir, filename), 'wb')
        output.write(response.read())
        output.close()

class Coral(Rutgers):
    """Setup Indonesial flowthrough"""
    def __init__(self, **kwargs):
        super(Coral, self).__init__(**kwargs)

    def jd2filename(self,jd):
        if jd == None: jd = 731583 
        filename = "/coral_avg_%05i.nc" % (jd - 731365)
        return "%s/%s" % (self.datadir, filename)     
