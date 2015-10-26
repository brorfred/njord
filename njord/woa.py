import os, urlparse

import numpy as np
import pylab as pl
from netCDF4 import Dataset

import base
import gmtgrid

class Base(base.Grid):
    """Setup World Ocean Atlas instance"""
    def __init__(self, res="1deg", **kwargs):
        if ("25" in str(res)) |  ("25" in str(res)):
            self.res = "0.25"
        elif "5" in str(res):
            self.res = "5deg"
        else:
            self.res = "1.00"
        super(Base, self).__init__(**kwargs)
        
    def setup_grid(self):
        """Setup necessary variables for grid """
        gridfilename = os.path.join(self.datadir, self.filename("temp"))
        if not os.path.isfile(gridfilename):
            self.retrive_file(self.data_url("temp")+self.filename("temp"), gridfilename)
        with Dataset(gridfilename, 'r') as g:
            self.lat = g.variables['lat'][:]
            self.lon = g.variables['lon'][:]
            self.llon,self.llat = np.meshgrid(self.lon,self.lat)
            self.zlev = g.variables['depth'][:]

    def load(self, fldname): #,jd=731583,yr=0,mn=1,dy=1,hr=3):
        """ Load climatological fields"""
        filename = os.path.join(self.datadir, self.filename(fldname))
        if not os.path.isfile(filename):
            self.retrive_file(self.data_url(fldname)+self.filename(fldname), filename)
        with Dataset(filename) as nc:     
            fld =  np.squeeze(nc.variables[self.vdict[fldname][1]][:]).copy()
            fld[fld>9999] = np.nan
            setattr(self, fldname, fld)

    @property
    def vdict(self):
        """Add a dict with all variables defined"""
        if not hasattr(self, "_vdict"):
            self._vdict = {
                'salt': ['salinity',                    's_an', "s","",   "all"],
                'temp': ['temperature',                 't_an', "t","v2", "decav"],
                'no3':  ['nitrate',                     'n_an', "n","",   "all"],
                'po4':  ['phosphate',                   '',     "p","",   "all"],
                'sio':  ['silicate',                    '',     "s","",   "all"],
                'o2ct': ['oxygen',                      'o_an', "o","",   "all"],
                'o2rl': ['oxygen_saturation',           'O_an', "O","",   "all"],
                'aou':  ['apparent_oxygen_utilization', 'A_an', "A","",   "all"]
                }
        return self._vdict
    
    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            g = pycdf.CDF(self.gridfile)
            self.landmask = g.var('A_an')[0,0,:,:]
            self.landmask[self.landmask < 9999] = 1
            self.landmask[self.landmask > 9999] = 0
        return self._landmask
                
class Woa13(Base):
    """Setup World Ocean Atlas instance"""
    def __init__(self, res="1deg", **kwargs):
        self.ver = "woa13"  
        self.timespan = kwargs.get("timespan", "decav")
        super(Woa13, self).__init__(**kwargs)

    def data_url(self, fld):
        print self.cdfurl
        return self.cdfurl.format(self.vdict[fld][3], self.vdict[fld][0],
                                  self.vdict[fld][4], self.res)

    def data_param(self, fld):
        return self.cdfparam % (self.vdict[fld][0], self.timespan, self.res,
                                self.filename(fld))
            
    def filename(self, fld, mn=0):
        res = {"5deg":"5d", "1.00":"01","0.25":"04"}[self.res]
        return "%s_%s_%s%02i_%s%s.nc" % (self.ver, self.vdict[fld][4],
                                       self.vdict[fld][2], mn, res, self.vdict[fld][3]) 
