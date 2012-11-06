from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.io import netcdf_file, netcdf_variable

import base

class Mass(base.Grid):
    """ Baseclass for all rutgers projects """
    def __init__(self, **kwargs):
        if not hasattr(kwargs, 'i2'): kwargs['i2'] = 54
        super(Mass, self).__init__(**kwargs)
        self.add_mp()
        self.add_landmask()

    def add_landmask(self):
        self.landmask = np.isnan(self.depth)

    def add_utmxy(self):
        """Add catersian coordinates """
        ncH = netcdf_file(self.gridfile)
        self.utmx = ncH.variables['x'][:]
        self.utmy = ncH.variables['y'][:]

    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile, 'r')
        mat = np.genfromtxt(self.datadir + '/mass_grid.dat')
        self.llon = np.zeros((self.jmt,self.imt))
        self.llat = np.zeros((self.jmt,self.imt))
        for vec in mat:
            self.llon[vec[1]-1, vec[0]-1] = -vec[2]
            self.llat[vec[1]-1, vec[0]-1] = vec[3]
        self.depth = g.variables['depth'][:].copy()
        self.depth[self.depth < 0.01] = np.nan

        #self.Cs_r = g.variables['Cs_r'][:]
        dayvec = g.variables['time'][:].astype(np.float).copy()
        self.jdvec = np.array(pl.datestr2num('2000-12-31') + dayvec)

    def load(self,fldname, **kwargs):
        """ Load velocity fields for a given day"""
        if fldname == "uv":
            self.load('u',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.load('v',jd=jd, yr=yr, mn=mn, dy=dy, hr=hr)
            self.uv = np.sqrt(self.u[:,1:,:]**2 + self.v[:,:,1:]**2)
            return
        self._timeparams(**kwargs)
        tpos = np.nonzero(self.jdvec <= self.jd)[0].max()
        vc = {'uvel': ['u',    'ecom.cdf',      9.155553e-05,  0],
              'vvel': ['v',    'ecom.cdf',      9.155553e-05,  0],
              'wvel': ['v',    'ecom.cdf',      6.103702e-08,  0],
              'temp': ['temp', 'ecom.cdf',      0.0005340739, 12.5],
              'salt': ['salt', 'ecom.cdf',      0.0006103702, 20],
              'chlo': ['chl',  'bem_water.cdf', 0.001525925,  50],
              'newp': ['np',   'bem_water.cdf', 0.001525925,  50],
              'netp': ['pp',   'bem_water.cdf', 0.001525925,  50],
              'tpoc': ['tpoc', 'bem_water.cdf', 0.01525925,  500],
              }
        try:
            nc = netcdf_file(self.datadir + vc[fldname][1])
        except KeyError:
            raise KeyError, "%s is not included" % fldname

        fld = nc.variables[vc[fldname][0]][tpos,:,
                                           self.j1:self.j2,
                                           self.i1:self.i2]
        fld = (fld * vc[fldname][2] + vc[fldname][3]).astype(np.float32)
        fld[:,self.landmask] = np.nan
        self.__dict__[fldname] = fld
        return tpos


        #self.ssh =  np.squeeze(nc.variables['zeta'][:])
        #self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
        #         self.Cs_r[:,np.newaxis,np.newaxis])
