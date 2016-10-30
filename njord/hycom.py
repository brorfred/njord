import os
from urlparse import urljoin
from datetime import datetime as dtm


import numpy as np
import pylab as pl
from scipy.io import netcdf_file
import netCDF4

import base
import gmtgrid

class GLBa0_08(base.Grid):
    """Setup Ecco 0.25 deg fields."""
    def __init__(self, **kwargs):
        self.pardict = {'uvel':'UVEL',
                        'vvel':'VVEL',
                        'wvel':'WVEL',
                        'salt':'SALT',
                        'temp':'THETA'}
        super(GLBa0_08, self).__init__(**kwargs)
        self.k1 = getattr(self, "k1", 0)
        self.k2 = getattr(self, "k2", self.kmt)
        self.add_mp()

    def setup_grid(self):
        """Setup necessary variables for grid """
        if not os.path.isfile(self.gridfile):
            self.download(self.gridfile, 'vvel')
        g = netcdf_file(self.gridfile, 'r')
        self.llat = g.variables['Latitude'][:]
        self.gmt = gmtgrid.Shift(g.variables['Longitude'][1649,:].copy())
        self.llon = self.gmt.field(g.variables['Longitude'][:].copy())
        self.llon[self.llon>180]  = self.llon[self.llon>180]-360
        self.llon[self.llon<-180] = self.llon[self.llon<-180]+360
        #self.llon[self.llon>180]  = np.nan
        #self.llon[self.llon<-180] = np.nan

        
        #self.zlev = g.variables['DEPTH_T'][:]
        #mat = np.squeeze(g.variables['VVEL'][:].copy())
        #mat[mat==0] = np.nan
        #self.depth = self.gmt.field(np.nanmax((
        #    self.zlev[:,np.newaxis,np.newaxis] * (mat*0+1)),axis=0))

    def load(self, fldname, **kwargs):
        """ Load velocity fields for a given day
        if fldname == "uv":
            self.load('uvel', **kwargs)
            self.load('vvel', **kwargs)
            self.uv = np.sqrt(self.uvel**2 + self.vvel**2)/2
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

        self.__dict__[fldname] = fld[self.k1:self.k2, self.j1:self.j2,
                                     self.i1:self.i2]
    """

    def download(self,filename, fieldname):
        """Download a missing file from source website
        print "Downloading file from server. This might take several minutes."
        url = urljoin(self.dataurl, "%s.nc" % self.pardict[fieldname],
                      os.path.basename(filename))
        self.retrive_file(url, filename)
    """ 
       

    @property
    def landmask(self):
        """ Generate a landmask for t-positions"""
        if not hasattr(self, '_landmask'):
            nc = netCDF4.Dataset(self.gridfile)
            
            self._landmask = self.gmt.field(nc.variables["u"][0,0,:,:].mask)


            
        return self._landmask


class GLBu0_08(base.Grid):
    """Setup Ecco 0.25 deg fields."""
    def __init__(self, **kwargs):
        self.pardict = {'uvel':'u', 'vvel':'v', 'salt':'SALT', 'temp':'THETA'}
        super(GLBu0_08, self).__init__(**kwargs)
        self.k1 = getattr(self, "k1", 0)
        self.k2 = getattr(self, "k2", self.kmt)
        self.add_mp()

    def setup_grid(self):
        """Setup necessary variables for grid """
        if not os.path.isfile(self.gridfile):
            self.download(self.gridfile, 'vvel')
        g = netcdf_file(self.gridfile, 'r')
        self.latvec = g.variables['lat'][:]
        self.gmt = gmtgrid.Shift(g.variables['lon'][:].copy())
        self.lonvec = self.gmt.lonvec
        self.llon,self.llat = np.meshgrid(self.lonvec,self.latvec)

        #self.zlev = g.variables['DEPTH_T'][:]
        #mat = np.squeeze(g.variables['VVEL'][:].copy())
        #mat[mat==0] = np.nan
        #self.depth = self.gmt.field(np.nanmax((
        #    self.zlev[:,np.newaxis,np.newaxis] * (mat*0+1)),axis=0))

    @property
    def landmask(self):
        """ Generate a landmask for t-positions"""
        if not hasattr(self, '_landmask'):
            nc = netCDF4.Dataset(self.gridfile)
            self._landmask = self.gmt.field(nc.variables["water_u"][0,0,:,:].mask)
        return self._landmask
