import os, os.path
import subprocess as sbp
#from datetime import datetime as dtm
import urllib
#import urllib2
#import re

import numpy as np
import pylab as pl
import netCDF4

import base

ZIPCMD = "pbzip2"

class AVHRR(base.Grid):

    def __init__(self, **kwargs):
        super(AVHRR, self).__init__(**kwargs)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.only_npz = kwargs.get('only_npz', False)

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        if not os.path.isfile(self.gridfile):
            print gridfile
            self.download(self.defaultjd)
            self._try_to_unzip(gridfile)
        gr = netCDF4.Dataset(self.gridfile)
        self.lonvec = gr.variables["lon"][:].copy()
        self.latvec = gr.variables["lat"][:].copy()
        self.lonvec[self.lonvec>180] = self.lonvec[self.lonvec>180]-360
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
        self.llon = np.roll(self.llon,720)
        
    def gen_filename(self, jd):
        """Generate filename"""
        self._timeparams(jd=jd)
        filestub = "%i/AVHRR/avhrr-only-v2.%i%02i%02i.nc"
        return filestub % (self.yr, self.yr, self.mn, self.dy)
            
    def refresh(self, jd1=None, jd2=None):
        """ Read a L3 mapped file and add field to current instance"""
        jd1 = self.jdmin if jd1 is None else jd1
        jd2 = int(pl.date2num(dtm.now())) - 1  if jd2 is None else jd2
        for jd in np.arange(jd1, jd2+1):
            filename = os.path.join(self.datadir, self.generate_filename(jd))
            print " --- %s --- " % pl.num2date(jd).strftime('%Y-%m-%d')
            print "Checking %s" % filename + '.bz2'
            if not os.path.isfile(filename + '.bz2'):
                try:
                    self.load(jd=jd, verbose=True)
                except IOError:
                    print "Downloading failed. Trying to remove old files."
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
                print "\n"
            else:
                print "found"
             
    def load(self, fld="sst", **kwargs):
        """Load the satellite field associated with a given time."""
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f" % self.jd)
        self.filename = os.path.join(self.datadir,self.gen_filename(self.jd))
        self.vprint( "Filename is %s" % (self.filename))
        """
        zipfname = self.filename + '.bz2'
        if not (os.path.isfile(self.filename) | os.path.isfile(zipfname)):
            print "File missing, downloading from GSFC"
            if not self.download(self.jd):
              print "File not on server, bailing"
              if hasattr(self, 'sst'):
                  del self.sst
              return False
        try:
            zipped = self._try_to_unzip(self.filename)
        except:
            return False
        """
        nc = netCDF4.Dataset(self.filename)
        nc.set_auto_mask("False")
        field = np.roll(np.squeeze(nc.variables["sst"]), 720)
        field[field<-100] = np.nan
        self.sst = (np.squeeze(field) / 10.)[self.j1:self.j2, self.i1:self.i2]
        #self._try_to_zip(zipped)
        return True

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
            gr = netcdf_file(self.gridfile)
            setattr(self, '_landmask', np.squeeze(gr.variables['mask'][:]==2))
        return self._landmask
                
    def download(self, jd):
        """Download a missing file from NODC"""
        self._timeparams(jd=jd)
        url = "%s/%04i/%03i/" % (self.dataurl, self.yr, self.yd)
        filename = self.generate_filename(jd) + ".bz2"
        self.vprint("Downloading %s" % url+filename)
        return self.retrive_file(url+filename, self.datadir + filename)

    def _try_to_unzip(self, filename):
        """Unzip file if exists and and is a valid bzip2 file"""
        if not os.path.isfile(filename):
            zipfname = filename + ".bz2"
            self.vprint( "Trying to uncompress file %s" % zipfname)
            err = sbp.call([ZIPCMD, "-d", zipfname])
            if err == 1:
                print "Decompression of " + zipfname + " failed."
                print "Trying to download again"
                self.download(self.jd)
                err = sbp.call([ZIPCMD, "-d", zipfname])
                if err == 1:
                    raise IOError, "Download file failed."
            return True
        else:
            return False

    def _try_to_zip(self, zipped):
        if zipped:
            self.vprint( "Compressing file")
            err = sbp.call([ZIPCMD, self.filename])
            if err ==1 :
                raise IOError( "Compression of " + self.filename + " failed.")




