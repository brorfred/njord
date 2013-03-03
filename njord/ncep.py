import os, urllib2
from datetime import datetime as dtm

import numpy as np
import pylab as pl
from scipy.io import netcdf_file

import base
import gmtgrid

class Daily(base.Grid):
    """ NCAR/NCEP daily reanalysis fields"""
    def __init__(self, **kwargs):
        super(Daily, self).__init__(**kwargs)
        self.add_mp()
        self.pardict = {'airp':'slp',}

    def setup_grid(self):
        """Setup necessary variables for grid """
        gc = netcdf_file(self.datadir + '/land.nc')
        self.lat = gc.variables['lat'][:]
        self.gmt = gmtgrid.Shift(gc.variables['lon'][:].copy())
        self.lon = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.landmask = self.gmt.field(
            gc.variables['land'][:].copy()).astype(np.bool)

    def load(self, fldname, **kwargs):
        """ Load NCEP reanalysis fields for a given day"""
        self._timeparams(**kwargs)
        filename = '%s/%s.%04i.nc' % (self.datadir,self.pardict[fldname],self.yr)
        if not os.path.exists(filename): self.download(filename)
        nc  = netcdf_file(filename)        
        fobj = nc.variables[self.pardict[fldname]]
        fld = self.gmt.field(fobj.data * fobj.scale_factor + fobj.add_offset)
        self.__dict__[fldname] = fld[self.yd-1,
                                     self.j1:self.j2, self.i1:self.i2]

    def download(self, filename):
        """Download a missing file from server"""
        print "Downloading file from server. This might take several minutes."
        remote_filename = os.path.basename(filename) 
        url = "%s/%s/%s" % (self.remote_domain, self.remote_path, remote_filename)
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(filename, 'wb')
        output.write(response.read())
        output.close()

    def add_landmask(self):
        """Load and add landmask"""
        pass
