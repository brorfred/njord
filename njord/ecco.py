import os, urllib2
from datetime import datetime as dtm

import numpy as np
import pylab as pl
from scipy.io import netcdf_file

import base
import gmtgrid

class Glob_025_ll(base.Grid):
    """Setup Ecco 0.25 deg fields."""
    def __init__(self, **kwargs):
        super(Glob_025_ll, self).__init__(**kwargs)
        self.add_mp()
        self.pardict = {'uvel':'UVEL',
                        'vvel':'VVEL',
                        'wvel':'WVEL',
                        'salt':'SALT',
                        'temp':'THETA'}

    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile, 'r')
        self.lat = g.variables['LATITUDE_T'][:]
        self.gmt = gmtgrid.Shift(g.variables['LONGITUDE_T'][:].copy())
        self.lon = self.gmt.lonvec 
        self.zlev = g.variables['DEPTH_T'][:]
        mat = np.squeeze(g.variables['VVEL'][:].copy())
        mat[mat==0] = np.nan
        self.depth = self.gmt.field(np.nanmax((
            self.zlev[:,np.newaxis,np.newaxis] * (mat*0+1)),axis=0))
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        
    def load(self, fldname, **kwargs):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u', **kwargs)
            self.load('v', **kwargs)
            self.uv = np.sqrt(self.u**2 + self.v**2)/2
            return
        self._timeparams(**kwargs)
        self.jd = np.round((self.jd+1)/3)*3-1
        self._jd_to_dtm()
        filename = ("%s/%s.1440x720x50.%04i%02i%02i.nc" %
                    (self.datadir, self.pardict[fldname],
                     self.yr,self.mn,self.dy))
        if not os.path.exists(filename):
            self.download(filename, fldname)
        nc = netcdf_file(filename)     
        fld =  self.gmt.field(np.squeeze(
            nc.variables[self.pardict[fldname]][:]).copy())
        fld[fld==0]  = np.nan
        fld[fld<-1e5] = np.nan

        self.__dict__[fldname] = fld[...,self.j1:self.j2, self.i1:self.i2]

    def download(self,filename, fieldname):
        """Download a missing file from GSFC's website"""
        print "Downloading file from server. This might take several minutes."
        uri = "http://ecco2.jpl.nasa.gov"
        folder = "data1/cube/cube92/lat_lon/quart_90S_90N/%s.nc" % self.pardict[fieldname]
        url = "%s/%s/%s" % (uri, folder, os.path.basename(filename))
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(filename, 'wb')
        output.write(response.read())
        output.close()

    @property
    def landmask(self):
        """ Generate a landmask for t-positions"""
        if not hasattr(self, '_landmask'):
            self._landmask = np.isnan(self.depth)
        return self._landmask
