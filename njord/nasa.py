import os, os.path
import glob
import subprocess as sbp
from datetime import datetime as dtm
import re
from distutils import spawn
import warnings
from urllib.parse import urlparse

from bs4 import BeautifulSoup                
import requests
import numpy as np
import njord.utils.mpl_dates as pl
from netCDF4 import Dataset

from njord import base
#import base
from njord.utils import yrday

class Base(base.Grid):

    def __init__(self, **kwargs):
        super(Base, self).__init__(**kwargs)
        self.datadir = '%s/%s%s/' % (self.datadir, self.fp, self.res)
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)
        if not hasattr(self, "maxjd"):
            self.maxjd = int(pl.date2num(dtm.now())) - 3

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        if self.res == "9km":
            i0t,imt,j0t,jmt = (0000 ,4320, 0, 2160)
        elif self.res == "4km":
            i0t,imt,j0t,jmt = (0000 ,8640, 0, 4320)
        elif self.res == "2km":
            i0t,imt,j0t,jmt = (0000 ,17280, 0, 8640)
        elif self.res == "1deg":
            i0t,imt,j0t,jmt = (0000 ,360, 0, 180)
        incr  = 360.0/imt
        jR    = np.arange(j0t, jmt)
        iR    = np.arange(i0t, imt)
        self.latvec = (  90 - jR*incr - incr/2)
        self.lonvec = (-180 + iR*incr + incr/2)
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
        self.imt = imt
        self.jmt = jmt
        for pstr, pval in zip(["i1","i2","j1","j2"], [i0t,imt,j0t,jmt]):
            setattr(self, pstr, getattr(self, pstr, pval))

    def generate_filename(self, fld='chl', fldtype="DAY", **kwargs):
        """Generate filename"""
        stamp = "%s%s.L3m_%s_%s_%skm.nc"
        if len(kwargs):
            self._timeparams(**kwargs)
        ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                 pl.date2num(dtm(self.yr,  1,  1))) + 1
        if fldtype.lower() == "mc":
            if not hasattr(self, "mc_datedict"):
                self.add_mnclim()
            datestr = self.mc_datedict[self.mn]
        elif "mo" in fldtype.lower():
            self.add_filepreflist(fldtype="mo")
            datestr = self.mo_fileprefs[self.yr*100 + self.mn]
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
            try:
                datestr = self._retrieve_datestamps(self.cu_url_9km)[-1]
            except (requests.ConnectTimeout,requests.ConnectionError):
                flist = glob.glob(os.path.join(self.datadir,
                    stamp % (self.fp, "*", fldtype, self.vc[fld], self.res[0])))
                if len(flist) == 0:
                    raise IOError("Can't download CU file")
                else:
                    return os.path.basename(flist[-1])
        else:
            raise TypeError(f"File type {fldtype} not implemented")
        return stamp % (self.fp, datestr, fldtype, self.vc[fld], self.res[0])

    def _retrieve_datestamps(self, path, ext=".nc"):
        """Retrive datestamps for all files in dir on GSFC server"""
        url = "https://%s/%s" % (self.gsfc_host, path)
        datelist = []
        bs = BeautifulSoup(
            requests.get(url, timeout=1).content, "html.parser")
        for td in bs.find_all("td"):
            if td.a is not None:
                datelist.append(
                    (urlparse(td.a.attrs["href"]).path.split("/")[-1][1:15]))
        if len(datelist) == 0:
            raise RuntimeError(f"No files found at {url}")
        return datelist

    def refresh(self, fld, fldtype="DAY", jd1=None, jd2=None, delall=False):
        """ Read a L3 mapped file and add field to current instance"""
        jd1 = self.minjd if jd1 is None else jd1
        jd2 = self.maxjd if jd2 is None else jd2 + 1
        for jd in np.arange(jd1, jd2):
            print(" --- %s --- " % pl.num2date(jd).strftime('%Y-%m-%d'))
            filename = os.path.join(
                self.datadir, self.generate_filename(fld,fldtype, jd=jd))
            if delall:
                for fn in glob.glob(filename + "*"):
                    print("Deleted %s" % fn)
                    os.remove(fn)
            print("Looking for file")
            if os.path.isfile(filename):
                print("Found")
            else:
                print("File doesn't exist.")
                self.download(filename)
             
    def load(self, fld="chl", fldtype="DAY", nan=np.nan, **kwargs):
        """Load the satellite field associated with a given time."""
        if fld in ["fqy",]:
            pass
        self._timeparams(**kwargs)
        self.vprint( "load jd=%f and field=%s" % (self.jd, fld))
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
        self.filename = os.path.join(
            self.datadir, self.generate_filename(fld, fldtype))
        self.vprint( "datadir:        %s" % self.datadir)
        self.vprint( "filename:       %s" % os.path.basename(self.filename))
        if not os.path.isfile(self.filename):
            status = self.download(self.filename)
            print(status)
        setattr(self, fld, self._l3read_nc4(fld, nan))

    @property
    def landmask(self):
        """Return a landmask"""
        if not hasattr(self,'_landmask'):
            try:
                fld = self.get_field('par', fldtype='CU', nan="nan")
            except KeyError:
                fld = self.get_field('chl', fldtype='CU', nan="nan")
            self._landmask = np.isnan(fld)
        return self._landmask
    
    def add_mnclim(self):
        """Add a list with file name dates for Monthly climatologies"""
        """
        datelist = self._retrieve_datestamps(getattr(self, f"mc_url_{self.res}"))
        self.mc_datedict = {}
        for dstr in datelist[1:]:
            mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
            self.mc_datedict[mn] = dstr
        """
        tvec = np.array(
            self._retrieve_datestamps(getattr(self, f"mc_url_{self.res}")))
        yrvec = np.array([int(tstr[0:4]) for tstr in tvec]) 
        ydvec = np.array([int(tstr[4:7]) for tstr in tvec])
        mnvec,_ = yrday.mndy(yrvec, ydvec) 
        self.mc_datedict = {}
        for mn in range(1,13):
            mnpart = tvec[mnvec==mn]
            maxpos = np.argmax([t[7:14] for t in mnpart])
            self.mc_datedict[mn] = mnpart[maxpos]

    def add_filepreflist(self, fldtype="MC"):
        """Add a list with file name dates for Monthly climatologies"""
        if hasattr(self, "%s_fileprefs" % fldtype):
            return
        predict = {}
        if "mc" in fldtype.lower():
            tvec = np.array(
                self._retrieve_datestamps(getattr(self, f"mc_url_{self.res}")))
            yrvec = np.array([int(tstr[0:4]) for tstr in tvec]) 
            ydvec = np.array([int(tstr[4:7]) for tstr in tvec])
            mnvec,_ = yrday.mndy(yrvec, ydvec) 
            self.mc_datedict = {}
            for mn in range(1,13):
                mnpart = tvec[mnvec==mn]
                maxpos = np.argmax([t[7:14] for t in mnpart])
                self.mc_datedict[mn] = mnpart[maxpos]
        elif "mo" in fldtype.lower():
            datelist = self._retrieve_datestamps(self.mo_url_9km)
            for dstr in datelist[1:]:
                mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
                yr = int(dstr[:4])
                predict[yr*100 + mn] = dstr
        setattr(self, "%s_fileprefs" % fldtype, predict)
        
    def fileurl(self, filename):
        return "%s%s" % (self.dataurl, os.path.basename(filename))
        
    def download(self, filename):
        """Download a missing file from GSFC's website"""
        self.retrive_file(self.fileurl(filename), local_filename=filename)
        if not os.path.isfile(filename):
            print("------------- Download failed! -------------""")

    def _l3read_nc4(self, fld, nan=np.nan):
        self.vprint( "Reading netCDF4 file")
        self.vprint(self.filename)
        with Dataset(self.filename) as nc:
            nc.set_auto_mask(True)
            nc.set_auto_scale(True)
        
            var = nc.variables[list(nc.variables.keys())[0]]
            raw = var[self.ijslice]
            if not hasattr(raw, "mask"):
                return self.llon * np.nan
            field = raw.data
            field[raw.mask] = np.nan  
            
        #start_jd    = pl.datestr2num(nc.time_coverage_start)
        #end_jd      = pl.datestr2num(nc.time_coverage_end)
        #self.jd     = ((start_jd + end_jd)/2)
        #self.date   = pl.num2date(self.jd)
        return field

    @property
    def fieldlist(self):
        return self.vc.keys()

class MODIS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "A"
        super(MODIS, self).__init__(**kwargs)

    @property
    def maxjd(self):
        return int(pl.date2num(dtm.now())) - 5

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return    {'k49':'KD490_Kd_490',     'adg':'IOP_adg_443_giop',
                   'chl':'CHL_chlor_a',      'poc':'POC_poc',
                   'bbp':'IOP_bbp_443_giop', 'sst':'SST_sst',
                   'par':'PAR_par', 'flh':'FLH_nflh', 'ipar':'FLH_ipar',
                   'eup':'ZLEE_Zeu_lee',     'pic':'PIC_pic'}

Aqua = MODIS
        
class SeaWiFS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res
        if res=="4km":
            raise TypeError("SeaWiFS has only 9km resolution")
        self.fp = "S"
        super(SeaWiFS, self).__init__(**kwargs)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return    {'k49':'KD490_Kd_490',
                   'chl':'CHL_chlor_a',
                   'poc':'POC_poc',
                   'bbp':'IOP_bbp_443_giop',
                   'par':'PAR_par',
                   'eup' :'ZLEE_Zeu_lee',
                   }


class OCTS(Base):
    def __init__(self, res="9km", **kwargs):
        if res=="4km":
            raise TypeError("SeaWiFS has only 9km resolution")
        self.res = res        
        self.fp = "O"
        super(OCTS, self).__init__(**kwargs)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return    {'k49':'KD490_Kd_490',
                   'chl':'CHL_chlor_a',
                   'poc':'POC_poc',
                   'bbp':'IOP_bbp_443_giop',
                   'par':'PAR_par',
                   }

        
class CZCS(Base):
    def __init__(self, res="9km", **kwargs):
        if res=="4km":
            raise TypeError("SeaWiFS has only 9km resolution")
        self.res = res
        self.fp = "C"
        super(CZCS, self).__init__(**kwargs)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return {'k49':'KD490_Kd_490', 'chl':'CHL_chlor_a'}

    
class VIIRS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "V"
        super(VIIRS, self).__init__(**kwargs)

    @property
    def vc(self):
        """Add a dict with filename variable components"""
        return    {'k49':'SNPP_KD490_Kd_490',
                   'chl':'SNPP_CHL_chlor_a',
                   'poc':'SNPP_POC_poc',
                   'bbp':'SNPP_IOP_bbp_443_giop',
                   'sst':'SNPP_SST_sst',
                   'par':'SNPP_PAR_par',
                   }


            
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
            print(jd,pl.num2date(jd).year,mn)
        return MCmin
