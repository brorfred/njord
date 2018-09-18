

import os 
import six

if six.PY2:
    raise ImportError("This version is not working wih python2 yet")



__all__ = ["ghrsst", "mimoc", "smos", "globcolour", "mldclim", 
           "avhrr", "gmtgrid", "nasa", "strang", "gompom", "fvcom", 
           "bem", "hycom", "ncep", "topaz", "bgccsm", "jpl", 
           "cesm", "landmask", "winds", "oscar", "woa", "ecco", "mati", 
           "poseidon", "worldclim", "gebco", "rutgers"]


from .cci    import OceanColor
from .ghrsst import L4_K10 as SST
