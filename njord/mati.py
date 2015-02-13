import os
from datetime import datetime as dtm
import urllib

import numpy as np
import pylab as pl

import numpy as np

from pyhdf.SD import SD,SDC
import base

class Cal(base.Grid):

    def __init__(self, **kwargs):
        super(Cal, self).__init__(**kwargs)
        self.add_mp()

    def setup_grid(self):
        """Setup necessary variables for grid """

        if not os.path.isfile(self.datadir + self.gridfile):
            urllib.urlretrieve(self.dataurl + self.gridfile,
                               self.datadir + self.gridfile)
        g = SD(self.datadir + self.gridfile, SDC.READ)
        self.llat = g.select('Latitude')[:]
        self.llon = g.select('Longitude')[:]
        
    def load(self, fldname, **kwargs):
        """ Load Cali Current fields for a given day"""
        self._timeparams(**kwargs)
        
        if fldname == 'chl':
            filename = "/C%04i%03i_chl_mapped.hdf" % (self.yr, self.yd)
            #ncfieldname = 'chl_%04i_%03i' % (yr,yd)
            def scale(PV): return 10**(PV*0.015-2)
        elif fldname == 'sst':
            filename = "/M%04i%03i_sst_mapped.hdf" % (self.yr, self.yd)
            #ncfieldname = 'sst_%04i_%03i' % (yr,yd)            
            def scale(PV): return PV*0.15000001-3
        if not os.path.isfile(self.datadir + filename):
            print "Downloading " + filename
            self.download(fldname, self.jd)
            
        h = SD(self.datadir + filename,SDC.READ)        
        ncfieldname = h.datasets().keys()[0]
        fld =  h.select(ncfieldname)
        attr = fld.attributes()
        PV = fld[:].astype(np.float)
        PV[PV<0] = PV[PV<0]+256
        PV[PV==0]   = np.nan
        PV[PV==255] = np.nan
        setattr(self, fldname, scale(PV)[self.j1:self.j2, self.i1:self.i2])


    def download(self, fldname, jd):
        self._timeparams(jd=jd)
        if fldname == "chl":
            url = "%s/%04i/C%04i_chl_day/" % (self.dataurl, self.yr, self.yr)
            filename = "/C%04i%03i_chl_mapped.hdf" % (self.yr, self.yd)
        elif fldname == "sst":
            url = "%s/%04i/M%04i_sst_day/" % (self.dataurl, self.yr, self.yr)
            filename = "/M%04i%03i_sst_mapped.hdf" % (self.yr, self.yd)
        else:
            raise ValueError("Wrong fieldname")
        urllib.urlretrieve(url+filename, self.datadir + filename)

    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            filename = os.path.basename(self.landurl)
            if not os.path.isfile(self.datadir + filename):
                urllib.urlretrieve(self.dataurl + self.landurl,
                                   self.datadir + filename)
            h = SD(self.datadir + filename, SDC.READ)        
            ncfieldname = h.datasets().keys()[0]
            self._landmask =   h.select(ncfieldname)[:] == -1
        return self._landmask
