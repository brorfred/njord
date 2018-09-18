
import shutil
import os, os.path
import subprocess as sbp
from datetime import datetime as dtm
import urllib
#import urllib2
#import re

import numpy as np
import njord.utils.mpl_dates as pl
import netCDF4

from . import base

if shutil.which("pbzip2"):
    ZIPCMD = "pbzip2"
elif shutil.which("bzip2"):
    ZIPCMD = "bzip2"
else:
    raise OSError("No bzip2 decompressor installed")
    
class Base(base.Grid):

    def __init__(self, **kwargs):
        self.maxjd = int(pl.date2num(dtm.now())) - 7
        super(Base, self).__init__(**kwargs)
        self.fieldlist = ["temp",]

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        self.gridfile = self.get_filename(jd=self.defaultjd)
        if not os.path.isfile(self.gridfile):
            self.download(self.defaultjd)
            self._try_to_unzip(self.gridfile)
        gr = netCDF4.Dataset(self.gridfile)
        self.lonvec = gr.variables['lon'][:].copy()
        self.latvec = gr.variables['lat'][:].copy()
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    def refresh(self, jd1=None, jd2=None):
        """ Read a L3 mapped file and add field to current instance"""
        if type(jd1) is str:
            jd1 = pl.datestr2num(jd1)
        elif jd1 is None:
            jd1 = self.jdmin
        if type(jd2) is str:
            jd2 = pl.datestr2num(jd2)
        elif jd1 is None:
            jd2 = int(pl.date2num(dtm.now()))
        for jd in np.arange(jd1, jd2+1):
            filename = os.path.join(self.datadir, self.get_filename(jd=jd))
            print(" --- %s --- " % pl.num2date(jd).strftime('%Y-%m-%d'))
            print("Checking %s" % filename + '.bz2')
            if not os.path.isfile(filename + '.bz2'):
                try:
                    self.load(jd=jd, verbose=True)
                except IOError:
                    print("Downloading failed. Trying to remove old files.")
                    try:
                        os.remove(filename)
                    except:
                        pass
                    try:
                        os.remove(filename + ".bz2")
                    except:
                        pass
                    try:
                        self.load(jd=jd,verbose=True)
                    except:
                        print ("   ###   Warning! Failed to add %s   ###" %
                               os.path.basename(filename))
                print("\n")
            else:
                print("found")
             
  

    def timestat(self, jd1, jd2):
        """Create a time-mean of field"""
        sumarr = np.zeros((12,) + self.llat.shape) * np.nan
        cntarr = np.zeros((12,) + self.llat.shape) * np.nan
        for jd in np.arange(jd1,jd2+1):
            self.vprint(" === Processing field %i of %i ===" % (jd-jd1, jd2-jd1))
            if not self.load(jd=jd):
                continue
            sumarr[self.mn-1,:,:] = np.nansum([sumarr[self.mn-1,:,:], self.sst],axis=0) 
            cntarr[self.mn-1,:,:] = np.nansum([cntarr[self.mn-1,:,:], self.sst*0+1],axis=0) 
        self.meansst = sumarr / cntarr

    @property
    def landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'_landmask'):
            gr = netCDF4.Dataset(self.gridfile)
            setattr(self, '_landmask', np.squeeze(gr.variables['mask'][:] == 2))
        return self._landmask.data[self.ijslice]
                
    def download(self, jd):
        """Download a missing file from NODC"""
        self._timeparams(jd=jd)
        url = "%s/%04i/%03i/" % (self.dataurl, self.yr, self.yd)
        filename = self.get_filename(jd=jd) + ".bz2"
        return self.retrive_file(url, filename)

    def _try_to_unzip(self, filename):
        """Unzip file if exists and and is a valid bzip2 file"""
        if not os.path.isfile(filename):
            zipfname = filename + ".bz2"
            self.vprint( "Trying to uncompress file %s" % zipfname)
            err = sbp.call([ZIPCMD, "-d", zipfname])
            if err == 1:
                print("Decompression of " + zipfname + " failed.")
                print("Trying to download again")
                self.download(self.jd)
                err = sbp.call([ZIPCMD, "-d", zipfname])
                if err == 1:
                    raise IOError("Download file failed.")
            return True
        else:
            return False

    def _try_to_zip(self, zipped):
        if zipped:
            self.vprint( "Compressing file")
            err = sbp.call([ZIPCMD, self.filename])
            if err ==1 :
                raise IOError( "Compression of " + self.filename + " failed.")

class L4_K10(Base):
    """Class to work with GHRSST L4 data"""
    def __init__(self, **kwargs):
        """Initiate class instance"""
        super(L4_K10, self).__init__(**kwargs)

    def get_filename(self, jd):
        """Generate filename"""
        self._timeparams(jd=jd)
        filestub = "%04i%02i%02i-NAVO-L4HR1m-GLOB-v01-fv01_0-K10_SST.nc"
        return os.path.join(self.datadir, filestub % (self.yr,self.mn,self.dy))

    def load(self, fld="temp", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f" % self.jd)
        self.filename = self.get_filename(jd=self.jd)
        self.vprint( "Filename is %s" % (self.filename))
        zipfname = self.filename + '.bz2'
        if not (os.path.isfile(self.filename) | os.path.isfile(zipfname)):
            self.vprint("File missing, downloading from server")
            if not self.download(self.jd):
              print("File not on server, bailing")
              if hasattr(self, 'sst'):
                  del self.sst
              return False
        try:
            zipped = self._try_to_unzip(self.filename)
        except:
            return False
        nc = netCDF4.Dataset(self.filename)
        raw = nc.variables['analysed_sst'][:].astype(np.float)
        fld = raw.data - 273.15
        fld[raw.mask] = np.nan
        setattr(self, "temp", np.squeeze(fld)[self.j1:self.j2,self.i1:self.i2])
        #self._try_to_zip(zipped)
        return True

class L4CMC(Base):
    """Class to work with GHRSST L4 data"""
    def __init__(self, **kwargs):
        """Initiate class instance"""
        super(L4CMC, self).__init__(**kwargs)

    def get_filename(self, jd):
        """Generate filename"""
        self._timeparams(jd=jd)
        filestub = ("%04i%02i%02i120000-CMC-L4_GHRSST-SSTfnd-" +
                    "CMC0.2deg-GLOB-v02.0-fv02.0.nc")
        return os.path.join(
            self.datadir, filestub % (self.yr,self.mn,self.dy))

    def download(self, jd):
        """Download a missing file from NODC"""
        self._timeparams(jd=jd)
        url = "%s/%04i/%03i/" % (self.dataurl, self.yr, self.yd)
        try:
            self.retrive_file(url, self.get_filename(jd=jd))
        except OSError as err:
            print(url)
            print(self.get_filename(jd=jd))
            raise OSError(err, url)
        
    def load(self, fld="temp", **kwargs):
        """Load the satellite field associated with a given time."""
        self._model_timeparams(**kwargs)
        self.vprint( "load jd=%f" % self.jd)
        self.filename = self.get_filename(jd=self.jd)
        self.vprint( "Filename is %s" % (self.filename))
        if not os.path.isfile(self.filename):
            try:
                self.download(self.jd)
                if hasattr(self, 'sst'):
                    del self.sst
            except OSError:
                print("Download failed")
                return False
        nc = netCDF4.Dataset(self.filename)
        fld  = nc.variables['analysed_sst'][0,...] - 273.15
        mask = nc.variables['mask'][0,...]
        fld[mask!=1] = np.nan
        setattr(self, "temp", np.squeeze(fld)[self.ijslice])


class SEVIRIatl(Base):
    """Class to work with GHRSST L4 data"""
    def __init__(self, **kwargs):
        """Initiate class instance"""
        super().__init__(**kwargs)

    def get_filename(self, **time_kwargs):
        """Generate filename"""
        if len(time_kwargs)>0:
            self._timeparams(**time_kwargs)
        if self.jd:
            self.meteosat_number = 10
        filename = (f"{self.yr:04}{self.mn:02}{self.dy:02}" + 
                    f"{self.hr:02}{self.min:02}00" +
                     "-OSISAF-L3C_GHRSST-SSTsubskin-SEVIRI_SST-ssteqc_" +
                    f"meteosat{self.meteosat_number:02}_" +
                    f"{self.yr:04}{self.mn:02}{self.dy:02}_" +
                    f"{self.hr:02}{self.min:02}00" +
                     "-v02.0-fv01.0.nc")
        return os.path.join(self.datadir, filename)

    def download(self, jd):
        """Download a missing file from NODC"""
        self._timeparams(jd=jd)
        url = "%s/%04i/%03i/" % (self.dataurl, self.yr, self.yd)
        try:
            self.retrive_file(url, self.get_filename(jd=jd))
        except OSError as err:
            print(url)
            print(self.get_filename(jd=jd))
            raise OSError(err, url)
        
    def load(self, fld="temp", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f" % self.jd)
        self.filename = self.get_filename(jd=self.jd)
        self.vprint( "Filename is %s" % (self.filename))
        if not os.path.isfile(self.filename):
            try:
                self.download(self.jd)
                if hasattr(self, 'sst'):
                    del self.sst
            except OSError:
                print("Download failed")
                return False
        nc = netCDF4.Dataset(self.filename)
        raw  = nc.variables['sea_surface_temperature'][0,...] - 273.15
        fld = raw.data
        fld[raw.mask] = np.nan
        setattr(self, "temp", np.squeeze(fld)[self.ijslice])
