
import os

import numpy as np
import pandas as pd
import netCDF4

from njord import base

class FayMcKinley(base.Grid):
    """Retrieve and load ecoregion definitions based on Fay and McKinley
    
    Description:
    https://doi.org/10.5194/essd-6-273-2014
    
    DOI to data:
    https://doi.org/10.1594/PANGAEA.828650
    
    Data url:
    https://doi.pangaea.de/10013/epic.42948.d019
    http://epic.awi.de/34786/19/Time_Varying_Biomes.nc
    """
    def __init__(self, filename=None):
        super().__init__()
            
    @property
    def regions(self):
        if not hasattr(self, "__regions"):
            self.load()
        return self.__regions

    @regions.setter
    def regions(self, regions):
        self.__regions = regions
    
    def setup_grid(self):
        if not os.path.isfile(self.filename):
            self.download()
        nc = netCDF4.Dataset(self.filename)
        self.latvec = nc.variables["lat"][:].data
        self.lonvec = nc.variables["lon"][:].data
        self.llon,self.llat = np.meshgrid(self.lonvec,self.latvec)
        self.jmt,self.imt = self.llon.shape

    def download(self):
        """Download the file"""
        filename = os.path.basename(self.filename)
        url = f"{self.dataurl}{filename}"
        self.retrive_file(url, local_filename=self.filename)
        
    def load(self, yr=None):
        """Load Biome array"""
        nc = netCDF4.Dataset(os.path.join(self.datadir, self.filename))
        fld = nc.variables["MeanBiomes"][:].data.T
        self.regions = fld
