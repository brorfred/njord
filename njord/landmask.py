import os
from distutils import spawn
import subprocess as sbp


import numpy as np

import requests
from netCDF4 import Dataset

import base

if spawn.find_executable("pbzip2"):
    ZIPCMD = "pbzip2"
elif spawn.find_executable("bzip2"):
    ZIPCMD = "bzip2"
else:
    raise IOError, "Couldn't find a bzip2 executable. Needed to unzip files"

class Navo1km(base.Grid):
    """Use the NAVO 1km landmask"""
    def __init__(self, **kwargs):
        super(Navo1km, self).__init__(**kwargs)

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        local_filename = os.path.join(self.datadir, self.gridfile)
        if not os.path.isfile(local_filename):
            print "Downloading gridfile"
            self.download()
        nc = Dataset(local_filename)
        self.latvec = nc.variables['lat'][:]
        self.lonvec = nc.variables['lon'][:]
        self.llon,self.llat = np.meshgrid(self.lonvec,self.latvec)
        self.landmask = nc.variables['dst'][:]        
                
    def download(self):
        """Download landmask file"""
        local_filename = os.path.join(self.datadir, self.gridfile + ".bz2")
        payload = {"m":"documents", "f":self.gridfile + ".bz2"}
        r = requests.get(url, params=payload, stream=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024): 
                if chunk:
                    f.write(chunk)
                    f.flush()
        try:
            err = sbp.call([ZIPCMD, "-d", local_filename])
        except OSError as det:
            raise OSError, "Error when unzipping %s: %s" % (zipfname,det)
