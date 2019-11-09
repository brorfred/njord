import os, os.path
import glob
import subprocess as sbp
from datetime import datetime as dtm
from distutils import spawn
import warnings
from pathlib import Path

import numpy as np

import netCDF4

from njord import nasa, base
#import nasa, base
import njord.utils.mpl_dates as pl

class Base(nasa.Base):
    def __init__(self, res="4km", **kwargs):
        self.res = res        
        self.fp = "ESACCI"
        self.maxjd = int(pl.date2num(dtm.now())) - 3
        super().__init__(**kwargs)
        self.timetypelist = ["5day", "8day", "daily", "monthly"]
        
    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        self.datadir = '%s/%s/' % (self.datadir, self.res)
        self.gridfile = os.path.join(self.datadir, self.generate_filename())
        self.vprint( "Grid datadir:  %s" % self.datadir)
        self.vprint( "Grid filename: %s" % os.path.basename(self.gridfile))
        if not os.path.isfile(self.gridfile):
            self.download()
        nc = netCDF4.Dataset(self.gridfile)
        self.lonvec = nc.variables["lon"][:].data
        self.latvec = nc.variables["lat"][:].data
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
        
    def download(self, fldname='chl', timetype="8D", **kwargs):
        """Download a missing file from the DAAC"""
        self._timeparams(**kwargs)
        filename = os.path.join(
            self.datadir, self.generate_filename(fldname, timetype))
        longfld = self.vc[fldname][0]
        tdict = self._get_timetype_dict(timetype)
        url = f"{self.dataurl}/{tdict['timetype']}/{longfld}/{self.yr}/"
        self.retrive_file(url, local_filename=filename)
        
    def load(self, fldname="chl", timetype="8day", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f and field=%s" % (self.jd, fldname))
        self.filename = os.path.join(
            self.datadir, self.generate_filename(fldname, timetype))
        self.vprint( "datadir:  %s" % self.datadir)
        self.vprint( "filename: %s" % os.path.basename(self.filename))
        if not self.isncfile(self.filename):
            self.download(fldname=fldname, timetype=timetype)
        nc = netCDF4.Dataset(self.filename)        
        raw = nc.variables[
            self.vc[fldname][0]][(slice(0,1),) + self.ijslice]
        field = raw.data
        field[raw.mask] = np.nan  
        setattr(self, fldname, np.squeeze(field))

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return    {#'k49': ['KD490_Kd_490',  'km', 0,      100  ],
                   'chl':  ['chlor_a',       'km', 0.001,  200  ],
                   #'ipar': ['FLH_ipar',     'km', 0,    10000  ],
                   #'eup': ['ZLEE_Zeu_lee',  'km', 0,    10000  ],
                   #'pic': ['PIC_pic',       'km', 0,        0.1],
                   }


class OceanColor(Base):
    def __init__(self, res="4km", **kwargs):
        super().__init__(**kwargs)
        self.timetypelist = ["5day", "8day", "daily", "monthly"]
        
    def _get_timetype_dict(self, raw_timetype):
        if "5" in raw_timetype:
            timetype = "5day"
            ttype = "5D_DAILY"
            jd5dvec = np.arange(0,365,5) + pl.date2num(dtm(self.yr, 1, 1))
            dt = pl.num2date(jd5dvec[self.jd>=jd5dvec].max())
            datestr = dt.strftime("%Y%m%d")
        elif "8" in raw_timetype:
            timetype = "8day"
            ttype = "8D_DAILY"
            jd8dvec = np.arange(0,365,8) + pl.date2num(dtm(self.yr, 1, 1))
            dt = pl.num2date(jd8dvec[self.jd>=jd8dvec].max())
            datestr = dt.strftime("%Y%m%d")
        elif "m" in raw_timetype.lower():
            timetype = "monthly"
            ttype = "1M_MONTHLY"
            datestr = f"{self.yr}{self.mn:02}"
        else:
            timetype = "daily"
            ttype = "1D_DAILY"
            datestr = f"{self.yr}{self.mn:02}{self.dy:02}"
        return {"timetype":timetype, "ttype":ttype, "datestr":datestr}
        
    def download(self, fldname='chl', timetype="8D", **kwargs):
        """Download a missing file from the DAAC"""
        self._timeparams(**kwargs)
        filename = os.path.join(
            self.datadir, self.generate_filename(fldname, timetype))
        longfld = self.vc[fldname][0]
        tdict = self._get_timetype_dict(timetype)
        url = f"{self.dataurl}/{tdict['timetype']}/{longfld}/{self.yr}/"
        self.retrive_file(url, local_filename=filename)

    def generate_filename(self, fldname='chl', timetype="8D", **kwargs):
        """Generate filename"""
        self._timeparams(**kwargs)
        ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                 pl.date2num(dtm(self.yr,  1,  1))) + 1
        self.filestamp= "%s-OC-L3S-%s-MERGED-%s_%s_GEO_PML_OCx-%s-%s.nc"
        tdict = self._get_timetype_dict(timetype)
        return(self.filestamp % (self.fp.upper(), self.vc[fldname][0].upper(),
                                 tdict["ttype"], self.res, tdict["datestr"],
                                 self.dataversion))

class LocalPML(Base):
    def __init__(self, res="4km", **kwargs):
        super().__init__(**kwargs)
        
    def _get_timetype_dict(self, raw_timetype):
        if "5" in raw_timetype:
            timetype = "5day"
            ttype = "5D_DAILY"
            jd5dvec = np.arange(0,365,5) + pl.date2num(dtm(self.yr, 1, 1))
            dt = pl.num2date(jd5dvec[self.jd>=jd5dvec].max())
            datestr = dt.strftime("%Y%m%d")
        elif "8" in raw_timetype:
            timetype = "8day"
            ttype = "8D_DAILY"
            jd8dvec = np.arange(0,365,8) + pl.date2num(dtm(self.yr, 1, 1))
            dt = pl.num2date(jd8dvec[self.jd>=jd8dvec].max())
            datestr = dt.strftime("%Y%m%d")
        elif "m" in raw_timetype.lower():
            timetype = "monthly"
            ttype = "1M_MONTHLY"
            datestr = f"{self.yr}{self.mn:02}"
        else:
            timetype = "daily"
            ttype = "1D_DAILY"
            datestr = f"{self.yr}{self.mn:02}{self.dy:02}"
        return {"timetype":timetype, "ttype":ttype, "datestr":datestr}

    def generate_filename(self, fldname='chl', timetype="8D", **kwargs):
        """Generate filename"""
        if ("m" in timetype.lower()) & ("c" in timetype.lower()):
            return ("ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_" +
                    "PML_OCx_QAA-%02i-fv3.1.nc" % self.mn)
        else:
            self._timeparams(**kwargs)
            ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                     pl.date2num(dtm(self.yr,  1,  1))) + 1
            self.filestamp= "%s-OC-L3S-%s-MERGED-%s_%s_GEO_PML_OCx-%s-%s.nc"
            tdict = self._get_timetype_dict(timetype)
            return(self.filestamp % (self.fp.upper(),
                                     self.vc[fldname][0].upper(),
                                     tdict["ttype"],
                                     self.res,
                                     tdict["datestr"],
                                     self.dataversion))

    
    def download(self, fldname='chl', timetype="8D", **kwargs):
        """Download a missing file from the DAAC"""
        self._timeparams(**kwargs)
        if ("m" in timetype.lower()) & ("c" in timetype.lower()):
            globaldir = self.mc_datadir
        else:
            longfld = self.vc[fldname][0]
            tdict = self._get_timetype_dict(timetype)
            globaldir = os.path.join(self.global_dir, tdict['timetype'],
                                     longfld, str(self.yr))
        filename = self.generate_filename(fldname, timetype)
        globalfile = os.path.join(globaldir, filename)
        os.symlink(globalfile, os.path.join(self.datadir, filename))
        


class LocalPMLSST(base.Grid):

    def __init__(self, **kwargs):
        if not hasattr(self, "maxjd"):
            self.maxjd = int(pl.date2num(dtm.now())) - 7
        super().__init__(**kwargs)
        self.fieldlist = ["sst", "temp"]

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        gr = netCDF4.Dataset(self.generate_filename(self.defaultjd))
        self.lonvec = gr.variables["lon"][:].copy()
        self.latvec = gr.variables["lat"][:].copy()
        self.lonvec[self.lonvec>180] = self.lonvec[self.lonvec>180]-360
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
        #self.llon = np.roll(self.llon,720)
        
    def generate_filename(self, jd):
        """Generate filename"""
        self._timeparams(jd=jd)
        datadir = os.path.join(
            self.datadir, f"{self.yr:04}/{self.mn:02}/{self.dy:02}/")
        filename = (f"{datadir}/{self.yr:04}{self.mn:02}{self.dy:02}120000-" +
                    "ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_" +
                    "CDR2.1-v02.0-fv01.0.nc")
        return os.path.join(datadir, filename)
             
    def load(self, fld="sst", **kwargs):
        """Load the satellite field associated with a given time."""
        if fld == "temp":
            self.load(fld="sst", **kwargs)
            self.temp = self.sst
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f" % self.jd)
        self.filename = os.path.join(self.generate_filename(self.jd))
        self.vprint( "Filename is %s" % (self.filename))
        nc = netCDF4.Dataset(self.filename)
        raw = nc.variables["analysed_sst"][:]
        fld = raw.data
        fld[raw.mask] = np.nan
        #fld = np.roll(np.squeeze(fld), 720)
        self.sst = np.squeeze(fld)[self.j1:self.j2, self.i1:self.i2]
        return None

    @property
    def landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'_landmask'):
            gr = netcdf_file(self.gridfile)
            setattr(self, '_landmask', np.squeeze(gr.variables['mask'][:]==2))
        return self._landmask
    
