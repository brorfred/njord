import os, urllib2, urlparse
from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.io import netcdf_file, netcdf_variable
from pydap.client import open_url

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
            self.load('u',**kwargs)
            self.load('v',**kwargs)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        
        if self.opendap:
            tpos = int(self.jd) - 714800
            k1   = kwargs.get("k1", getattr(self, "k1", self.klev)) 
            k2   = kwargs.get("k2", getattr(self, "k2", k1+1))
            dapH = open_url(self.dapurl)
            fld  = dapH[fldname][tpos,k1:k2,self.j1:self.j2,self.i1:self.i2] 
        else:
            filename = self.jd2filename(self.jd)
            if not os.path.isfile(filename):
                print "File missing"
                url = urlparse.urljoin(self.dataurl, os.path.basename(filename))
                self.retrive_file(url, filename)
            with netcdf_file(filename) as nc:
                fld =  nc.variables[fldname][:].copy()
                self.ssh =  np.squeeze(nc.variables['zeta'][:])
                self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
                              self.Cs_r[:,np.newaxis,np.newaxis])
        
        fld[fld>9999] = np.nan
        setattr(self, fldname, np.squeeze(fld))

                
class NWA(Rutgers):
    """Setup North West Atlantic"""
    def __init__(self, **kwargs):
        super(NWA, self).__init__(**kwargs)

    def jd2filename(self, jd):
        filename = "nwa_avg_%05i.nc" % (jd - 714782)
        return os.path.join(self.datadir, filename)    

class Coral(Rutgers):
    """Setup Indonesial flowthrough"""
    def __init__(self, **kwargs):
        super(Coral, self).__init__(**kwargs)

    def jd2filename(self,jd):
        if jd == None: jd = 731583 
        filename = "/coral_avg_%05i.nc" % (jd - 731365)
        return "%s/%s" % (self.datadir, filename)     
