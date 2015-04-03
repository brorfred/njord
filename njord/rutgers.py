import os, urllib2, urlparse
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

    @property
    def landmask(self):
        """ Return attribute """
        if not hasattr(self, "_landmask"):
            with netcdf_file(self.gridfile) as gH:
                self._landmask = gH.variables['mask_rho'][:]
        return self._landmask
            
    def add_utmxy(self):
        """Add catersian coordinates """
        g = netcdf_file(self.gridfile)
        self.utmx = g.variables['x_rho'][:]
        self.utmy = g.variables['y_rho'][:]

    def setup_grid(self):
        """Setup necessary variables for grid """
        if not os.path.isfile(self.gridfile):
            url = urlparse.urljoin(self.gridurl, os.path.basename(self.gridfile))
            self.retrive_file(url, self.gridfile)
        with netcdf_file(self.gridfile) as gH:
            self.llat  = gH.variables['lat_rho'][:]
            self.llon  = gH.variables['lon_rho'][:]-360
            self.depth = gH.variables['h'][:]
            self.Cs_r  = gH.variables['Cs_r'][:]

    def load(self,fldname, **kwargs):
        """ Load velocity fields for a given day"""
        self._timeparams(**kwargs)
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        
        filename = self.jd2filename(self.jd)
        if not os.path.isfile(filename):
            url = urlparse.urljoin(self.dataurl, os.path.basename(filename))
            print url
            print filename
            self.retrive_file(url, filename)
        with netcdf_file(filename) as nc:
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

    def jd2filename(self, jd):
        filename = "nwa_avg_%05i.nc" % (jd - 714782)
        return os.path.join(self.datadir, filename)    

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
