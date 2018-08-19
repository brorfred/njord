import os
import glob
from datetime import datetime as dtm

import numpy as np
import pylab as pl
import netCDF4

import base

class CSF3A(base.Grid):
    def __init__(self, **kwargs):
        super(CSF3A, self).__init__(**kwargs)

    def setup_grid(self):
        gc = netCDF4.Dataset(self.gridfile)
        self.latvec = gc.variables['lat'][:]
        self.lonvec = gc.variables['lon'][:]
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    @property
    def filelist(self):
        if not hasattr(self, "_filelist"):
            self._filelist = glob.glob(os.path.join(self.datadir, "*.nc"))
        return self._filelist

    @property
    def fulltvec(self):
        mask = (self.filejds[:,0]>=self.minjd) & (self.filejds[:,0]<=self.maxjd)
        return np.sort(self.filejds[:,0][mask])
        
    @property
    def filejds(self):
        if not hasattr(self, "_filejds"):
            jdlist = []
            for fn in self.filelist:
                jd1 = pl.datestr2num(os.path.basename(fn)[19:27])
                jd2 = pl.datestr2num(os.path.basename(fn)[35:43])
                jdlist.append([jd1,jd2,jd2-jd1])
            self._filejds = np.array(jdlist)
        return self._filejds

    @property
    def minjd(self):
        return self.filejds[:,0][self.filejds[:,2] <= 10].min()

    @property
    def maxjd(self):
        return self.filejds[:,1][self.filejds[:,2] <= 10].max()
    

    def filename_10d(self, jd):
        fpos = np.nonzero((self.filejds[:,0] <= jd) &
                          (self.filejds[:,2] <= 10))[0].max()
        return self.filelist[fpos]
            
    def load(self,fldname="salt", **kwargs):
        """ Load SMOS fields for a given day"""
        self._timeparams(**kwargs)

        nc = netCDF4.Dataset(self.filename_10d(self.jd))
        nc.set_auto_mask(False)
        fld = nc.variables["Mean_Sea_Surface_Salinity"][self.j1:self.j2,
                                                        self.i1:self.i2]
        fld[fld<0] = np.nan
        setattr(self, fldname, fld)

    @property
    def fieldlist(self):
        return ["salt",]

    def download(self, filename):
        """Download a missing file from source website"""
        print "Downloading file from server. This might take several minutes."
        url = urljoin(self.dataurl, os.path.basename(filename))
        self.retrive_file(self.dataurl, filename + ".gz")
        err = sbp.call([ZIPCMD, "-d", filename + ".gz"])


    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            self.load("uvel")
            self._landmask = np.isnan(self.uvel)
        return self._landmask
        
