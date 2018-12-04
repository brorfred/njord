import numpy as np
from osgeo import gdal, ogr    

class Shape2Array(object):

    def __init__(self, shapefile, layername=None,
                 attrfield=None, pixel_size=2.5, njord_inst=None, glob=False):
        self.orig  = ogr.Open(shapefile)
        self.layername = layername
        self.pixel_size = pixel_size
        self.njord_inst = njord_inst
        self.glob = glob
        self.get_layer()
        self.setup_grid()
        self.attrfield = (self.fieldnames[0]
                          if attrfield is None else attrfield)
                
    def setup_grid(self):
        nj = self.njord_inst
        if nj is not None:
            self.ymin = nj.llat.min()
            self.xmin = nj.llon.min()
            self.ymax = nj.llat.max()
            self.xmax = nj.llon.max()
            self.imt =  nj.imt * 2
            self.jmt =  nj.jmt * 2
            self.lat = np.linspace(self.ymin, self.ymax, self.jmt)
            self.lon = np.linspace(self.xmin, self.xmax, self.imt)
            self.llon,self.llat = np.meshgrid(self.lon, self.lat)
        else:
            if self.glob:
                (self.xmin, self.xmax,
                 self.ymin, self.ymax) = (-180,180,-90,90)
            else:
                (self.xmin, self.xmax,
                 self.ymin, self.ymax) = self.layer.GetExtent()
            self.imt = int((self.xmax-self.xmin) / self.pixel_size)
            self.jmt = int((self.ymax-self.ymin) / self.pixel_size)
            self.lat = np.linspace(self.ymin, self.ymax, self.jmt)
            self.lon = np.linspace(self.xmin, self.xmax, self.imt)
            self.llon,self.llat = np.meshgrid(self.lon, self.lat)

    def get_layer(self):
        if self.orig.GetLayerCount() == 1:
            self.layer = self.orig.GetLayer(0)
            self.layername = self.layer.GetName()
        elif self.layername is not None:
            self.layer = self.orig.GetLayerByName(self.layername)
        else:
            nlayer = np.arange(orig.GetLayerCount())
            errortext = ("More than one layer, please set " +
                         "'layername'. Available layers: %s" %
                         "\n".join([orig.GetLayerByIndex(n)
                                    .GetName() for n in nlayer]))
            raise KeyError(errortext)

    def create_mem_layer(self):
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
        
    def create_raster(self):
        self.raster = gdal.GetDriverByName('MEM').Create(
                      'arr', self.imt, self.jmt, 3, gdal.GDT_Byte)
        self.raster.SetGeoTransform((self.xmin, self.pixel_size, 0, 
                                     self.ymax, 0, -self.pixel_size))
        band = self.raster.GetRasterBand(1)
        band.SetNoDataValue(-999)
        if self.source_srs:
            self.raster.SetProjection(self.source_srs.ExportToWkt())
        else:
            # Source has no projection (needs GDAL >= 1.7.0 to work)
            self.raster.SetProjection('LOCAL_CS["arbitrary"]')

    def rasterize(self):
        self.create_mem_layer()
        self.add_numerical_attribute_fields()
        self.create_raster()
        err = gdal.RasterizeLayer(self.raster, [1],
                                  self.source_layer,
                                  options=["ATTRIBUTE=%s" %
                                           self.attrfield])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)
        arr = self.raster.GetRasterBand(1).ReadAsArray()
        arr[arr == 0] = -998
        self.array = arr - 1
        return self.array[::-1,:]
    
    def list_fields(self):
        for n,feature in enumerate(self.fieldlist):
            print("% 4i %s" % (n,feature))
