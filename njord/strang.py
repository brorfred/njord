import os, os.path
from urlparse import urljoin
from collections import OrderedDict
from cStringIO import StringIO

import numpy as np
import requests
try:
    import h5py
    HAS_H5PY = True
except:
    HAS_H5PY = False

import base
from utils import yrday

class After2005(base.Grid):

    def __init__(self, **kwargs):
        super(After2005, self).__init__(**kwargs)
    
    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""

        if not os.path.isfile(self.latfile):
            self.retrive_file(urljoin(self.gridurl,"la11.dat"), self.latfile)
        if not os.path.isfile(self.lonfile):
            self.retrive_file(urljoin(self.gridurl,"lo11.dat"), self.lonfile)
        self.llat = np.loadtxt(self.latfile)
        self.llon = np.loadtxt(self.lonfile)

                
    def load(self, fld, mn=None, **kwargs):
        """Load the data field associated with a given time."""
        self._timeparams(**kwargs)
        par = self.vc[fld][-1]
        filename = self.get_filename(par)
        if not os.path.isfile(filename):
            self.download(par)
        setattr(self, fld, np.loadtxt(filename))

                
    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'landmask'):
            self.load('par','CU', nan="nan")
            self.landmask = np.isnan(self.par)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return  {'uvirr':['CIE UV irradiance',      'mW/m^2',      116],
                 'glirr':['Global irradiance',      'W/m^2',       117],
                 'diirr':['Direct irradiance',      'W/m^2',       118],
                 'durr': ['Sunshine duration',      'min',         119],
                 'par':  ['Photosynthetic photons', 'mumol/s/m^2', 120],
                }

    def get_filename(self, par, lev=0):
        parlist = [par, self.yr, self.mn, self.dy, self.hr, lev]
        filename = "_".join(["%04i" % v for v in parlist]) + ".asc"
        return os.path.join(self.datadir, filename)

    def _get_h5filename(self, fld):
        fn = "%s_%04i%02i%02i.h5" % (fld, self.yr, self.mn, self.dy) 
        return os.path.join(self.datadir, fn)

    def write_h5(self, fld, arr, lev=0):
        """Write field to hdf5 file."""
        with h5py.File(self._get_h5filename(fld, lev=lev), 'a') as h5f:
            if not fld in h5f:
                h5f.create_dataset(fld, (24,) + self.shape)
            h5f.get(fld)[self.hr,:,:] = arr

    def read_h5(self, fld, lev=0):
        """Write field to hdf5 file."""
        with h5py.File(self._get_h5filename(fld, lev=lev), 'r') as h5f:
            if fld in h5f:
                return h5f.get(fld)[self.hr,:,:]
            else:
                return False
            
    def isfield(self, fld):
        """Check if a data field is stored locally"""
        if HAS_H5PY & os.path.isfile(self._get_h5filename(fld)):
            pass
        else:
            return os.path.isfile(self.get_filename(par))
            
    def download(self, fld, lev=0):
        """Download missing data files"""
        print "Downloading datafiles. This will take a while."
        payload = OrderedDict([('par', self.vc[fld][-1]), ('y1', self.yr),
                               ('m1', self.mn), ('d1', self.dy),
                               ('h1', self.hr), ('lev', lev)])
        if HAS_H5PY:
            tempfile = os.path.join(tempdir, "strang.asc")
            for h in np.arange(24):
                payload['h1'] = h + 1
                arr = np.loadtxt(StringIO(self.retrive_file(self.dataurl,
                                                            params=payload)))
                self.write_h5(fld, arr, lev=lev)
        else:
            self.retrive_file(self.dataurl, self.get_filename(self.vc[par][-1]),
                              params=payload)

        
    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            self._landmask =  ~np.isnan(self.get_field('alt'))
        return self._landmask        
