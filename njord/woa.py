from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.spatial import KDTree, cKDTree
import matplotlib.cm as cm
from scipy.io import netcdf_file, netcdf_variable

import base
import gmtgrid

class Woa09(base.Grid):
    """Setup World Ocean Atlas instance"""
    def __init__(self, res="1deg", **kwargs):
        self.res = res
        super(Woa09, self).__init__(**kwargs)
        self.add_mp()
        self.add_vdict()
        
    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file("%s_%s.nc" % (self.gridfile,self.res), 'r')
        self.lat = g.variables['lat'][:]
        self.gmt = gmtgrid.Shift(g.variables['lon'][:].copy())
        self.lon = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        #self.depth = g.variables['depth'][:]

    def load(self,fldname ): #,jd=731583,yr=0,mn=1,dy=1,hr=3):
        """ Load climatological fields"""
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        filename = "/%s_annual_%s.nc" % (self.vdict[fldname][0],self.res)
        nc = netcdf_file("%s/%s" % (self.datadir, filename))     
        fld =  np.squeeze(nc.variables[self.vdict[fldname][1]][:]).copy()
        fld[fld>9999] = np.nan
        self.__dict__[fldname] = self.gmt.field(fld)


    def add_vdict(self):
        """Add a dict with all variables defined"""
        self.vdict = {
            'salt': ['salinity', 's_an'],
            'temp': ['temperature', 't_an'],
            'no3':  ['nitrate', 'n_an'],
            'po4':  ['phosphate', ''],
            'sio':  ['silicate', ''],
            'o2ct': ['dissolved_oxygen', 'o_an'],
            'o2rl': ['oxygen_saturation', 'O_an'],
            'aou':  ['apparent_oxygen_utilization', 'A_an']
            }

    def add_landmask(self):
        g = pycdf.CDF(self.gridfile)
        self.landmask = g.var('A_an')[0,0,:,:]
        self.landmask[self.landmask < 9999] = 1
        self.landmask[self.landmask > 9999] = 0

    def add_utmxy(self):
        g = pycdf.CDF(self.datadir + '/coral_grd.nc')
        self.utmx = g.var('x_rho')[:]
        self.utmy = g.var('y_rho')[:]
        

