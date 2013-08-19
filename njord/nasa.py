import os, os.path
import subprocess as sbp
from datetime import datetime as dtm
import urllib
import urllib2
import re

import numpy as np
import pylab as pl

from pyhdf.SD import SD, SDC
import base
import yrday

class MODIS(base.Grid):

    def __init__(self, res="9km", **kwargs):
        self.res = res
        super(MODIS, self).__init__(**kwargs)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.add_vc()
        self.datadir = self.datadir + '/A' + self.res + '/' 
        
    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        if self.res == "9km":
            self.i1,self.i2,self.j1,self.j2 = (0000 ,4320, 0, 2160)
        elif self.res == "4km":
            self.i1,self.i2,self.j1,self.j2 = (0000 ,8640, 0, 4320)
        incr  = 360.0/self.i2
        jR    = np.arange(self.j1, self.j2)
        iR    = np.arange(self.i1, self.i2)
        [x,y] = np.meshgrid(iR,jR)
        self.llat = (  90 - y*incr - incr/2)
        self.llon = (-180 + x*incr + incr/2)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.imt = self.i2
        self.jmt = self.j2

    def load(self,fld,fldtype="DAY", nan="nan", **kwargs):
        """Load the satellite field associated with a given time."""
        if fld in ["fqy",]:
            pass
        self._timeparams(**kwargs)
        if fld == 'Dchl':
            if getattr(self, 'djd1', 0) != self.jd:
                self.load('chl',fldtype, jd=self.jd)    
                self.dmat1 = self.chl
            self.load('chl',fldtype, jd=self.jd+1)
            self.dmat2 = self.chl
            self.Dchl = self.dmat2 - self.dmat1
            self.djd1 = self.jd+1
            self.dmat1 = self.dmat2
            return
     
        ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                 pl.date2num(dtm(self.yr,  1,  1))) + 1    
        if fldtype == "MC":
            self.add_mnclim()
            datestr = self.mc_datedict[self.mn]
        elif fldtype == "DAY":
            datestr = "%i%03i" % (self.yr, self.yd)
        elif fldtype == "8D":
            yd1 = np.arange(1,365,8)
            yd2 = np.arange(8,370,8)
            yd2[-1] = ydmax
            pos = np.nonzero(self.yd >= yd1)[0].max()
            datestr = ("%i%03i%i%03i" % 
                       (self.yr, yd1[pos], self.yr, yd2[pos]))
        elif fldtype == "CU":
            datestr = "20021852011365"
        else:
            print "Field type not included"
        self.filename = ("A%s.L3m_%s_%s_%s%s" % (datestr, fldtype,
                                            self.vc[fld][0],
                                            self.res[0],
                                            self.vc[fld][1]))
        self._l3read(fld,nan=nan)
        
    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if not hasattr(self,'landmask'):
            self.load('par','CU', nan="nan")
            self.landmask = np.isnan(self.par)
        
    def add_mnclim(self):
        """Add a list with file name dates for Monthly climatologies"""
        url = "http://%s/%s" % (self.gsfc_host, self.a_mc_url_9km)
        datelist = []
        for line in urllib.urlopen(url):
            if "cgi/getfile" in line:
                datelist.append(re.findall(r'getfile/A([^E]+).L3m',
                                           line)[0][:14])
        self.mc_datedict = {}
        for dstr in datelist[1:]:
            mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
            self.mc_datedict[mn] = dstr
        
    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['KD490_Kd_490','km',0,100],
                   'chl': ['CHL_chlor_a','km',0.001,200],
                   'poc': ['POC_poc','km',0,10000],
                   'bbp': ['GIOP01_bbp_443_giop','km',0],
                   'sst': ['SST','',-2,45],
                   'par': ['PAR_par','km',-5,200],
                   'flh': ['FLH_nflh','km',0],
                   'ipar': ['FLH_ipar','km',0],
                   'eup': ['KDLEE_Zeu_lee','km',0],
                   'eup_lee': ['KDLEE_Zeu_lee','km',0],
                   'eup_mor': ['KDMOREL_Zeu_morel','km',0]
                   }

    def download(self, filename):
        """Download a missing file from GSFC's website"""
        uri = "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/"
        url = "%s%s.bz2" % (uri, os.path.basename(filename))
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(filename + ".bz2",'wb')
        output.write(response.read())
        output.close()

    def _l3read(self,fld='',nan="nan"):
        """ Read a L3 mapped file and add field to current instance"""
        self.filename = self.datadir + '/' + self.filename
        if not fld: fld = self.filename.rsplit('_')[-2]

        if os.path.isfile(self.filename + '.npz'):
            l3m_data,base,intercept,slope = self._l3read_npz(fld)
            zipped = False
        else:
            zipfname = self.filename + '.bz2'
            if not (os.path.isfile(self.filename) | 
                    os.path.isfile(zipfname)):
                print "File missing, downloading from GSFC"
                self.download(self.filename)
            zipped = self._try_to_unzip(zipfname)
            zipped = True
            l3m_data,base,intercept,slope = self._l3read_hdf()

        if base != -999:
            self.__dict__[fld] = base**((slope * l3m_data) + intercept)
        else:
            self.__dict__[fld] = ((slope * l3m_data) + intercept)
        mask = ((self.__dict__[fld] < self.minval) |
                (self.__dict__[fld] > self.maxval))
        if nan=="nan":
            self.__dict__[fld][mask] = np.nan
        self._try_to_zip(zipped)
        if ((not os.path.isfile(self.filename + '.npz')) &
            ((self.llat.shape) == (self.jmt,self.imt))):
            self.add_ij()
            self._l3write_npz(l3m_data[~mask], self.imat[~mask],
                              self.jmat[~mask], base, intercept, slope)

    def _l3read_hdf(self, fieldname='l3m_data'):
        sd = SD(self.filename, SDC.READ)
        ds = sd.select(fieldname)
        field     = ds[self.j1:self.j2, self.i1:self.i2].copy()
        intercept = ds.attributes()['Intercept']
        slope     = ds.attributes()['Slope']
        try:
            nanval = ds.attributes()['Fill']
        except:
            nanval = 65535
        try:
            base  = ds.attributes()['Base']
        except KeyError:
            base = -999
        if 'Scaled Data Maximum' in sd.attributes().keys():
            self.maxval = sd.attributes()['Scaled Data Maximum']
            self.minval = sd.attributes()['Scaled Data Minimum']
        elif 'Suggested Image Scaling Maximum' in sd.attributes().keys():
            self.maxval = sd.attributes()['Suggested Image Scaling Maximum']
            self.minval = sd.attributes()['Suggested Image Scaling Minimum']
        else:
            self.minval = self.vc[fld][2]
            self.maxval = self.vc[fld][3]
        start_iso = (pl.date2num(dtm(
                        sd.attributes()['Period Start Year'],1,1)) + 
                        sd.attributes()['Period Start Day'] - 1)
        end_iso   = (pl.date2num(dtm(
                        sd.attributes()['Period End Year'],1,1)) + 
                        sd.attributes()['Period End Day'] - 1)
        self.jd   = ((start_iso+end_iso)/2)
        self.date = pl.num2date(self.jd)
        return field,base,intercept,slope

    def _l3read_npz(self, fld):
        """Read an npz file with field data as vectors""" 
        l3m_data = self.llat * 0 + self.vc[fld][2] - 1
        fH = np.load(self.filename + '.npz')
        dvec      = fH['dvec']
        ivec      = fH['ivec']
        jvec      = fH['jvec']
        base      = fH['base']
        intercept = fH['intercept']
        slope     = fH['slope']
        self.jd   = fH['jd']
        self.date = fH['date']
        mask = ((jvec>=self.j1) & (jvec<self.j2-1) &
                (ivec>=self.i1) & (ivec<self.i2-1))
        l3m_data[jvec[mask]-self.j1, ivec[mask]-self.i1] = dvec[mask]
        self.minval = self.vc[fld][2]
        self.maxval = self.vc[fld][3]
        return l3m_data,base,intercept,slope

    def _l3write_npz(self, dvec, ivec ,jvec, base, intercept, slope):
        """Create an npz file with field data as vectors""" 
        kwargs = {'dvec' : dvec,
                  'ivec' : ivec.astype(np.int16),
                  'jvec' : jvec.astype(np.int16),
                  'base' : base, 'intercept':intercept, 'slope':slope,
                  'jd'   : self.jd, 'date':self.date}
        np.savez_compressed(self.filename + '.npz', **kwargs)

    def _try_to_unzip(self, zipfname):
        """Unzip file if exists and and is a valid bzip2 file"""
        if not os.path.isfile(self.filename):
            err = sbp.call(["pbzip2", "-d", zipfname])
            if err == 1:
                print "Decompression of " + zipfname + " failed."
                print "Trying to download again"
                self.download(self.filename)
                err = sbp.call(["pbzip2", "-d", zipfname])
                if err == 1:
                    raise IOError, "Download file failed."
            return True
        else:
            return False

    def _try_to_zip(self, zipped):
        if zipped:
            err = sbp.call(["pbzip2", self.filename])
            if err ==1 :
                raise IOError( "Compression of " + self.filename + " failed.")

    def histmoller(self, fieldname, jd1, jd2, y1, y2,
                   mask=[], bins=100, logy=True):
        if len(mask)==0: mask = (self.lat !=-800)
        if logy:
            vlist = np.exp(np.linspace(np.log(y1), np.log(y2), bins+1))
        else:
            vlist = np.linspace(y1, y2, bins+1)
        hsmat = np.zeros((jd2-jd1+1,bins), dtype=np.int)
        tvec = np.arange(jd1,jd2+1)
        for n_t,jd in enumerate(tvec):
            print pl.num2date(jd)
            self.load(fieldname, jd=jd)
            field = self.__dict__[fieldname]
            field[~mask] = np.nan
            hsmat[n_t,:],_ = np.histogram(field[~np.isnan(field)], vlist)

        class hiS: pass
        hiS.tvec = tvec
        hiS.vlist = vlist
        hiS.hist = hsmat
        hiS.norm = (hsmat.astype(float) / hsmat.max(axis=0))
        self.__dict__[fieldname + "_hiS"] = hiS






def mission_min():
        i1 = 0
        i2 = None
        j1 = 0
        j2 = None
        t1 = pl.date2num(dtm(2003,1,1))
        t2 = pl.date2num(dtm(2010,12,31))

        ns = nasa()
        MCmin = np.zeros((12,ns.lat.shape[0],ns.lat.shape[1])) * 999

        for jd in np.arange(t1,t2):
            ns.load('chl',jd=jd)
            mn = pl.num2date(jd).month-1
            MCmin[mn,:,:] = np.nanmin(
                np.vstack( (MCmin[mn:+1,:,:],ns.chl[np.newaxis,:,:]) ),
                axis=0)
            print jd,pl.num2date(jd).year,mn
        return MCmin
