
import os
import glob

import numpy as  np
from geocube.api.core import make_geocube
import geopandas as pd
import xarray as xr

from njord import base

class Longhurst(base.Grid):
    """Retrieve and load ecoregion definitions based on Longhurst
    
    Description:
    https://doi.org/10.1093/plankt/17.6.1245
        
    Data url:
    http://oceandata.azti.es:8080
    /thredds/fileServer/MESMA/Longhurst_world_v4_2010.shp
    """
    def __init__(self, Dlat=0.1, Dlon=0.1, grid=None, filename=None, attrfield=None):
        if filename is None:
            datadir = os.path.dirname(__file__)
            self.filename = os.path.join(datadir, "data" ,"longhurst_v4_2010", "Longhurst_world_v4_2010.shp")
        else:
            self.filename = filename
        self.Dlon = Dlon
        self.Dlat = Dlat
        self._griddef = grid
        self.read_shp(self.filename)
        self._ds = self.shp2arr()
        super().__init__()
  
    def setup_grid(self):
        """Setup grid for the njord object"""
        self.latvec = self._ds.y.values[::-1]
        self.lonvec = self._ds.x.values[1:]
        self.llon,self.llat = np.meshgrid(self.lonvec,self.latvec)
        self.jmt,self.imt = self.llon.shape
        (self.xmin, self.xmax) = (self.lonvec.min(), self.lonvec.max())
        (self.ymin, self.ymax) = (self.latvec.min(), self.latvec.max())

    def read_shp(self, filename):
        """Read shape file into geopandas object"""
        self._geodf = pd.read_file(filename)
        self._geodf["longhurstID"] = self._geodf.index
        self.attrlist = list(self._geodf.keys())
        self.attrlist.remove("geometry")

    def shp2arr(self, geodf=None, Dlon=None, Dlat=None):
        "Rasterize geopandas dataframe to xarray dataset"
        geodf = self._geodf if geodf is None else geodf
        Dlon = self.Dlon if Dlon is None else Dlon
        Dlat = self.Dlat if Dlat is None else Dlat
        return make_geocube(self._geodf, resolution=(-Dlat, Dlon))

    def load(self, fldname="regions", *args, **kwargs):
        """Load Biome array"""
        if hasattr(self, "patch_array"):
            return None
        elif hasattr(self, "_region_file"):
            ds = xr.open_dataset(self._region_file)
            setattr(self, "_regions", ds.longhurstID.data)
        else:
            pass

    @property
    def regions(self):
        return self._ds.longhurstID.values[::-1,1:]

    @property
    def fieldnames(self):
        """Return a list of fields in the shapefile layer"""
        return list(self._geodf["ProvCode"])

    @property
    def region_names(self):
        """Get list of long names of all providences"""
        return list(self._geodf["ProvDescr"])

    @property
    def region_codes(self):
        """Get list of large-scale sections of providences"""
        return [nm.split("-")[0].rstrip() for nm in self.region_names]
