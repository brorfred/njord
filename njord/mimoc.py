"""Module to load and define fields from the different mixed layer
climatologies to the njord framework. The resulting field will be stored
as a attribute in the class instance. The depth  meters.  

"""
import os

import numpy as np
import pylab as pl
from scipy.io import netcdf_file

import gmtgrid
import base

class Brest(base.Grid):
    """ Mixed layer climatology from Montegut et al. See more information at
    http://www.locean-ipsl.upmc.fr/~clement/mld.html
    MLD can be defined in three different ways:
    Temperature (temp) - Depth where temp is  10m temp +- 0.2 deg C
    Density (dens) - Depth where dens is  10m dens + 0.03 kg m-3
    Variable criterion (var)
    
    """
    def __init__(self, type="dens", **kwargs):
        """Initialize with chosen type of MLD"""
        if type == "temp":
            self.mldFile = "mld_DT02_c1m_reg2.0.nc"
        elif type == "var":
            self.mldFile = "mld_DReqDTm02_c1m_reg2.0.nc"
        else:
            self.mldFile = "mld_DR003_c1m_reg2.0.nc"
        super(Brest, self).__init__(**kwargs)
        self.add_mp()

    def setup_grid(self):
        """Define lat and lon matrices for njord"""
        gc = netcdf_file(self.gridfile)
        self.lat = gc.variables['lat'][:].copy()
        self.gmt = gmtgrid.Shift(gc.variables['lon'][:].copy())
        self.lon = self.gmt.lonvec         
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        
    def load(self,mn=1):
        """ Load mixed layer climatology for a given day"""
        self.mn = mn
        fn = os.path.join(self.datadir, self.mldFile)
        if not os.path.isfile(fn):
            self.retrive_file(self.cdfurl+"/"+self.mldFile, local_filename=fn)
        nc = netcdf_file(fn)
        self.mld = self.gmt.field(nc.variables['mld']
                                  [self.mn-1, self.j1:self.j2, self.i1:self.i2])
        self.mld[self.mld< 0]=np.nan
        self.mld[self.mld>1e4]=np.nan

    def add_landmask(self):
        """Add landmask defined as 1=land 0=ocean"""
        gc = pycdf.CDF(self.datadir + "/" + self.gridfile)
        u = nc.var('u')[:,:,:,:1081]
        msk = np.squeeze(nanmean(u,axis=0))
        msk = msk*0
        msk[np.isnan(msk)]=1
        self.landmask = gmtgrid.convert(msk,self.gr)

