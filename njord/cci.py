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
        longfld = self.vc[fldname][1]
        tdict = self._get_timetype_dict(timetype)
        url = f"{self.dataurl}/{tdict['timetype']}/{longfld}/{self.yr}/"
        self.retrive_file(url, local_filename=filename)
        
    def load(self, fldname="chl", timetype="daily", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f and field=%s" % (self.jd, fldname))
        self.filename = os.path.join(
            self.datadir, self.generate_filename(fldname, timetype))
        self.vprint( "datadir:  %s" % self.datadir)
        self.vprint( "filename: %s" % os.path.basename(self.filename))
        if not self.isncfile(self.filename):
            self.download(fldname=fldname, timetype=timetype)
        with netCDF4.Dataset(self.filename) as nc:
            raw = nc.variables[
                self.vc[fldname][0]][(slice(0,1),) + self.ijslice]
            field = raw.data
            field[raw.mask] = np.nan  
            setattr(self, fldname, np.squeeze(field))

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return {'chl':         ['chlor_a',     'chlor_a', 'OCx'],
                'MODISA_nobs': ['MODISA_nobs', 'chlor_a', 'OCx'],
                'MERIS_nobs':  ['MERIS_nobs',  'chlor_a', 'OCx'],
                'SeaWiFS_nobs':['SeaWiFS_nobs','chlor_a', 'OCx'],
                'VIIRS_nobs':  ['VIIRS_nobs',  'chlor_a', 'OCx'],
                'rrs412': ['Rrs_412', 'rrs', 'RRS'],
                'rrs443': ['Rrs_443', 'rrs', 'RRS'],
                'rrs490': ['Rrs_490', 'rrs', 'RRS'],
                'rrs510': ['Rrs_510', 'rrs', 'RRS'],
                'rrs555': ['Rrs_555', 'rrs', 'RRS'],
                'rrs670': ['Rrs_670', 'rrs', 'RRS'],
                #'ipar': ['FLH_ipar',     'km'],
                #'eup': ['ZLEE_Zeu_lee',  'km'],
                #'pic': ['PIC_pic',       'km'],
                }


class OceanColor(Base):
    def __init__(self, res="4km", ver = "3.1", **kwargs):
        self.data_version = f"{ver:.1f}" if type(ver) is not str else ver
        super().__init__(**kwargs)
        self.timetypelist = ["5day", "8day", "daily", "monthly"]

    def download(self, fldname='chl', timetype="8D", **kwargs):
        """Download a missing file from the DAAC"""
        self._timeparams(**kwargs)
        filename = os.path.join(
            self.datadir, self.generate_filename(fldname, timetype))
        longfld = self.vc[fldname][1]
        tdict = self._get_timetype_dict(timetype)
        url = f"{self.dataurl}/{tdict['timetype']}/{longfld}/{self.yr}/"
        self.retrive_file(url, local_filename=filename)

    def generate_filename(self, fldname='chl', timetype="8D", **kwargs):
        """Generate filename"""
        self._timeparams(**kwargs)
        ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                 pl.date2num(dtm(self.yr,  1,  1))) + 1
        self.filestamp= "%s-OC-L3S-%s-MERGED-%s_%s_GEO_PML_%s-%s-%s.nc"
        tdict = self._get_timetype_dict(timetype)
        return(self.filestamp % (self.fp.upper(), self.vc[fldname][1].upper(),
                                 tdict["ttype"], self.res, 
                                 self.vc[fldname][2], tdict["datestr"],
                                 self._data_version))

class LocalPML(Base):
    def __init__(self, res="4km", ver=3.1, **kwargs):
        self.data_version = f"{ver:.1f}" if type(ver) is not str else ver
        super().__init__(**kwargs)
        self.global_dir = self.global_dir.replace("3.1", self.data_version)
    
    def generate_filename(self, fldname='chl', timetype="8D", **kwargs):
        """Generate filename"""
        if ("m" in timetype.lower()) & ("c" in timetype.lower()):
            return ("ESACCI-OC-MAPPED-CLIMATOLOGY-1M_MONTHLY_4km_" +
                    f"PML_OCx_QAA-{self.mn:02}-fv{self.data_version}.nc")
        elif "month" in timetype.lower():
            return ("ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_" +
                   f"GEO_PML_OCx-{self.yr}{self.mn:02}-" +
                   f"fv{self.data_version}.nc")
        else:
            self._timeparams(**kwargs)
            self.filestamp= "%s-OC-L3S-%s-MERGED-%s_%s_GEO_PML_%s-%s-fv%s.nc"
            tdict = self._get_timetype_dict(timetype)
            return(self.filestamp % (self.fp.upper(),
                                     self.vc[fldname][1].upper(),
                                     tdict["ttype"],
                                     self.res,
                                     self.vc[fldname][2],
                                     tdict["datestr"],
                                     self.data_version))

    
    def download(self, fldname='chl', timetype="8D", **kwargs):
        """Download a missing file from the DAAC"""
        self._timeparams(**kwargs)
        if ("m" in timetype.lower()) & ("c" in timetype.lower()):
            globaldir = self.mc_datadir
        else:
            longfld = self.vc[fldname][1]
            tdict = self._get_timetype_dict(timetype)
            globaldir = os.path.join(self.global_dir, tdict['timetype'],
                                     longfld, str(self.yr))
            globaldir = globaldir.replace("v3.1",f"v{self.data_version}")
        filename = self.generate_filename(fldname, timetype)
        globalfile = os.path.join(globaldir, filename)
        os.symlink(globalfile, os.path.join(self.datadir, filename))
        
        #/data/datasets/CCI/v3.1-release/geographic/netcdf/monthly/chlor_a 


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
        with netCDF4.Dataset(self.filename) as nc:
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
    
