
from . import ogr2array

from .fay_mckinley import FayMcKinley

from .longhurst import Longhurst


def longhurst(Dlonlat=0.5, **kwargs):
    fn = "data/longhurst_v4_2010/Longhurst_world_v4_2010.shp"
    sh = ogr2array.Shape2Array(fn, pixel_size=Dlonlat, **kwargs)
    mat = sh.rasterize().astype(int)
    return {"regions":mat, "llon":sh.llon, "llat":sh.llat,
            "regioncodes":sh.get_fields()["ProvCode"],
            "regionnames":sh.get_fields()["ProvDescr"]}


