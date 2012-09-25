import os, os.path
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
        if self.res is "9km":
            self.i1,self.i2,self.j1,self.j2 = (0000 ,4320, 0, 2160)
        elif self.res is "4km":
            self.i1,self.i2,self.j1,self.j2 = (0000 ,8640, 0, 4320)
        incr  = 360.0/self.i2
        jR    = np.arange(self.j1, self.j2)
        iR    = np.arange(self.i1, self.i2)
        [x,y] = np.meshgrid(iR,jR)
        self.llat = (  90 - y*incr - incr/2)
        self.llon = (-180 + x*incr + incr/2)
        self.imt = self.i2
        self.jmt = self.j2

    def load(self,fld,fldtype="DAY",
             yr=2005, yd=300, mn=1, dy=0, jd=0, nan="nan"):
        """Load the satellite field associated with a given time."""
        if fld in ["fqy",]:
            pass
        
        if jd != 0:
            yr = pl.num2date(jd).year
            yd = jd - pl.date2num(dtm(yr,1,1)) + 1
            mn = pl.num2date(jd).month
        elif dy != 0:
            yd = (pl.date2num(dtm(yr,mn,dy)) -
                  pl.date2num(dtm(yr,1,1))) + 1
        elif yd < 1:
            yr = yr -1
            ydmax = (pl.date2num(dtm(yr, 12, 31)) -
                     pl.date2num(dtm(yr,  1,  1))) + 1    
            yd = ydmax + yd            
        ydmax = (pl.date2num(dtm(yr, 12, 31)) -
                 pl.date2num(dtm(yr,  1,  1))) + 1    
        if fldtype == "MC":
            self.add_mnclim()
            datestr = self.mc_datedict[mn]
        elif fldtype == "DAY":
            datestr = "%i%03i" % (yr, yd)
        elif fldtype == "8D":
            yd1 = np.arange(1,365,8)
            yd2 = np.arange(8,370,8)
            yd2[-1] = ydmax
            pos = np.nonzero(yd >= yd1)[0].max()
            datestr = "%i%03i%i%03i" % (yr, yd1[pos], yr, yd2[pos])
        elif fldtype == "CU":
            datestr = "20021852011365"
        else:
            print "Field type not included"
        filename = ("A%s.L3m_%s_%s_%s%s" % (datestr, fldtype,
                                            self.vc[fld][0],
                                            self.res[0],
                                            self.vc[fld][1]))
        self._l3read(filename,fld,nan=nan)
        
    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if hasattr(self,'landmask'):
            return
        self.load('par','CU')
        self.par = self.par * 0
        self.par[np.isnan(self.par)] = 1
        self.landmask = self.par.copy().astype(np.int8)
        del(self.par)
        
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
        self.vc = {'k49': ['KD490_Kd_490','km',0],
                   'chl': ['CHL_chlor_a','km',0],
                   'poc': ['POC_poc','km',0],
                   'bbp': ['GIOP01_bbp_443_giop','km',0],
                   'sst': ['SST','',0],
                   'par': ['PAR_par','km',0],
                   'flh': ['FLH_nflh','km',0],
                   'ipar': ['FLH_ipar','km',0],
                   'eup': ['KDLEE_Zeu_lee','km',0],
                   'eup_lee': ['KDLEE_Zeu_lee','km',0],
                   'eup_mor': ['KDMOREL_Zeu_morel','km',0]
                   }

    def download(self,filename):
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

    def _l2read(self ,filename,fld):
        i1,i2,j1,j2 = self.ijarea
        """
        try:
            os.system('pbzip2 -d ' + filename)
        except:
            print "Decompression of " + filename + " failed."
            print "File maybe not compressed."
            filename = filename[:-4]
        """
        fH = sd.select(fld)
        fld = fH[i1:i2,j1:j2].astype('float')
        fld[fld==fH.attributes()['bad_value_scaled']] = np.nan
        Intercept = fH.attributes()['intercept']
        Slope     = fH.attributes()['slope']
        try:
            Base  = fH.attributes()['base']
            fld = Base**((Slope*fld) + Intercept)
        except KeyError:
            fld = ((Slope*fld) + Intercept)
        """
        start_iso = (pl.date2num(datetime.datetime(
                    sd.attributes()['Period Start Year'],1,1)) +
                     sd.attributes()['Period Start Day'] - 1)
        end_iso   = (pl.date2num(datetime.datetime(
                    sd.attributes()['Period End Year'],1,1)) +
                     sd.attributes()['Period End Day'] - 1)
        self.date = pl.num2date((start_iso+end_iso)/2)
        self.iso  = ((start_iso+end_iso)/2)

        try:
            os.system('pbzip2 ' + filename)
        except:
            print "Compression of " + filename + " failed."
        """
        return fld

    def _l3read(self ,filename,fld='',nan="nan"):
        """ Read a L3 mapped file and add field to current instance"""
        filename = self.datadir + '/' + filename
        if not (os.path.isfile(filename) |
                os.path.isfile(filename +'.bz2')):
            print "File missing, downloading from GSFC"
            self.download(filename)
        if not fld: fld = filename.rsplit('_')[-2]

        zipped = False
        if os.path.isfile(filename + ".bz2"):
            try:
                os.system('pbzip2 -d ' + filename + ".bz2")
            except:
                raise IOerror( "Decompression of " + filename + " failed.")
            zipped = True

        sd = SD(filename, SDC.READ)
        l3m_data  = np.ma.masked_values(
            sd.select('l3m_data')[self.j1:self.j2, self.i1:self.i2],65535)
        Intercept = sd.attributes()['Intercept']
        Slope     = sd.attributes()['Slope']
        try:    
            Base  = sd.attributes()['Base']
            self.__dict__[fld] = Base**((Slope*l3m_data) + Intercept)
        except KeyError:
            self.__dict__[fld] = ((Slope*l3m_data) + Intercept)
        if nan=="nan":
            self.__dict__[fld][self.__dict__[fld] < self.vc[fld][2]] = np.nan
        else:
            self.__dict__[fld][self.__dict__[fld] < self.vc[fld][2]] = nan

        start_iso = (pl.date2num(dtm(
                    sd.attributes()['Period Start Year'],1,1)) + 
                     sd.attributes()['Period Start Day'] - 1)
        end_iso   = (pl.date2num(dtm(
                    sd.attributes()['Period End Year'],1,1)) + 
                     sd.attributes()['Period End Day'] - 1)
        self.date = pl.num2date((start_iso+end_iso)/2)
        self.jd   = ((start_iso+end_iso)/2)

        if zipped:
            try:
                os.system('pbzip2 ' + filename)
            except:
                raise IOerror( "Compression of " + filename + " failed.")

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
