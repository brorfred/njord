
import os
import glob

import numpy as  np
import pylab as  pl
import pandas as pd
from scipy.spatial import cKDTree

import numpy as np
from osgeo import gdal, ogr    

from njord import base


class Longhurst(base.Grid):
    """Retrieve and load ecoregion definitions based on Longhurst
    
    Description:
    https://doi.org/10.1093/plankt/17.6.1245
        
    Data url:
    http://oceandata.azti.es:8080
    /thredds/fileServer/MESMA/Longhurst_world_v4_2010.shp
    """
    def __init__(self, Dlatlon=0.1, filename=None, attrfield=None):
        self.Dlatlon = Dlatlon
        super().__init__()
        os.environ["SHAPE_RESTORE_SHX"] = "YES"
        self.orig  = ogr.Open(self.filename)
        self.layer = self.orig.GetLayer(0)
        self.layername = self.layer.GetName()
        layer_defn = self.layer.GetLayerDefn()
        print(self.filename)
        self.attrfield = (self.fieldnames[0]
                          if attrfield is None else attrfield)
                
    def setup_grid(self):
        if not os.path.isfile(self.filename):
            self.download()
        self.latvec = np.arange( -90, 90, self.Dlatlon)
        self.lonvec = np.arange(-180,180, self.Dlatlon)
        self.llon,self.llat = np.meshgrid(self.lonvec,self.latvec)
        self.jmt,self.imt = self.llon.shape
        (self.xmin, self.xmax) = (self.lonvec.min(), self.lonvec.max())
        (self.ymin, self.ymax) = (self.latvec.min(), self.latvec.max())


    def download(self):
        """Download the file"""
        for ext in ["shp", "dbf"]:
            filename = os.path.basename(self.filename)
            url = f"{self.dataurl}{filename}".replace("shp", ext)
            self.retrive_file(url, local_filename=self.filename.replace("shp", ext))
        
    def load(self, fldname="regions", jd=None):
        """Load Biome array"""
        if hasattr(self, "patch_array"):
            return None
        self._create_mem_layer()
        self.add_numerical_attribute_fields()
        self.create_raster()
        err = gdal.RasterizeLayer(
            self.raster, [1], self.source_layer,
            options=["ATTRIBUTE=%s" % self.attrfield])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)
        arr = self.raster.GetRasterBand(1).ReadAsArray()
        arr[arr == 0] = -998
        self.patch_array = arr - 1
        self._regions =  self.patch_array[::-1,:]

    @property
    def regions(self):
        if not hasattr(self, "_regions"):
            self.load()
        return self._regions

    def _create_mem_layer(self):
        self.source_ds = ogr.GetDriverByName("Memory").CopyDataSource(self.orig, "")
        self.source_layer = self.source_ds.GetLayer(0)
        self.source_srs = self.source_layer.GetSpatialRef()
    
    def add_numerical_attribute_fields(self):
        field_def = ogr.FieldDefn(self.attrfield, ogr.OFTReal)
        self.source_layer.CreateField(field_def)
        source_layer_def = self.source_layer.GetLayerDefn()
        field_index = source_layer_def.GetFieldIndex(self.attrfield)
        for idx,feature in enumerate(self.source_layer):
            feature.SetField(field_index, idx+1)
            self.source_layer.SetFeature(feature)

    def create_raster(self):
        self.raster = gdal.GetDriverByName('MEM').Create(
                      'arr', self.imt, self.jmt, 3, gdal.GDT_Byte)
        self.raster.SetGeoTransform((
            self.xmin, self.Dlatlon, 0, self.ymax, 0, -self.Dlatlon))
        band = self.raster.GetRasterBand(1)
        band.SetNoDataValue(-999)
        if self.source_srs:
            self.raster.SetProjection(self.source_srs.ExportToWkt())
        else:
            # Source has no projection (needs GDAL >= 1.7.0 to work)
            self.raster.SetProjection('LOCAL_CS["arbitrary"]')

    @property
    def fieldnames(self):
        layer_defn = self.layer.GetLayerDefn()
        return [layer_defn.GetFieldDefn(i).GetName()
                for i in range(layer_defn.GetFieldCount())]

    def get_fields(self):
        [setattr(self,fn,[]) for fn in self.fieldnames]
        self.layer.ResetReading()
        for ff in self.layer:
            for key in ff.keys():
                getattr(self,key).append(ff[key])
        return {key:getattr(self,key) for key in ff.keys()}
        
    @property
    def region_names(self):
        """Get list of long names of all providences"""
        return self.get_fields()["ProvDescr"]

    @property
    def region_codes(self):
        """Get list of large-scale sections of providences"""
        return [nm.split("-")[0].rstrip()
                    for nm in self.get_fields()["ProvDescr"]]
