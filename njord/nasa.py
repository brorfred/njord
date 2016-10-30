import os, os.path
import glob
import subprocess as sbp
from datetime import datetime as dtm
import urllib
import urllib2
import re
from distutils import spawn
import warnings

import numpy as np
import pylab as pl

from netCDF4 import Dataset

import base
from njord.utils import yrday

class Base(base.Grid):

    def __init__(self, **kwargs):
        super(Base, self).__init__(**kwargs)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.add_vc()
        self.datadir = '%s/%s%s/' % (self.datadir, self.fp, self.res) 
        self.only_npz = kwargs.get('only_npz', False)
        if not hasattr(self, "maxjd"):
            self.maxjd = int(pl.date2num(dtm.now())) - 3
        self.use_npz = setattr(self, "use_npz", False)

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        if self.res == "9km":
            i1,i2,j1,j2 = (0000 ,4320, 0, 2160)
        elif self.res == "4km":
            i1,i2,j1,j2 = (0000 ,8640, 0, 4320)
        elif self.res == "2km":
            i1,i2,j1,j2 = (0000 ,17280, 0, 8640)
        elif self.res == "1deg":
            i1,i2,j1,j2 = (0000 ,360, 0, 180)
        incr  = 360.0/i2
        jR    = np.arange(j1, j2)
        iR    = np.arange(i1, i2)
        [x,y] = np.meshgrid(iR,jR)
        self.llat = (  90 - y*incr - incr/2)
        self.llon = (-180 + x*incr + incr/2)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.imt = i2
        self.jmt = j2
        self.i1 = i1 if not hasattr(self, 'i1') else self.i1
        self.i2 = i2 if not hasattr(self, 'i2') else self.i2
        self.j1 = j1 if not hasattr(self, 'j1') else self.j1
        self.j2 = j2 if not hasattr(self, 'j2') else self.j2

    def generate_filename(self, fld='chl', fldtype="DAY", **kwargs):
        """Generate filename"""
        if len(kwargs):
            self._timeparams(**kwargs)
        ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                 pl.date2num(dtm(self.yr,  1,  1))) + 1
        if fldtype == "MC":
            self.add_mnclim()
            datestr = self.mc_datedict[self.mn]
        if "mo" in fldtype.lower():
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
            self.a_cu_url_9km = 'MODISA/Mapped/Cumulative/4km/chlor/'
            datestr = max(self._retrieve_datestamps(self.a_cu_url_9km))
        else:
            raise TypeError, "File average type not included"
        return("%s%s.L3m_%s_%s_%s%s.nc" % (self.fp, datestr, fldtype,
										self.vc[fld][0], self.res[0],
                                        self.vc[fld][1]))

    def _retrieve_datestamps(self, path, ext=".nc"):
        """Retrive datestamps for all files in dir on GSFC server"""
        url = "http://%s/%s" % (self.gsfc_host, path)
        datelist = []
        for line in urllib.urlopen(url):
            if ("cgi/getfile" in line) & (ext in line):
                datelist.append(re.findall(
                    r'getfile/%s([^E]+).L3m' % self.fp, line)[0][:14])
        return datelist


    def refresh(self, fld, fldtype="DAY", jd1=None, jd2=None, delall=False):
        """ Read a L3 mapped file and add field to current instance"""
        jd1 = pl.datestr2num('2003-01-01') if jd1 is None else jd1
        jd2 = int(pl.date2num(dtm.now())) - 1  if jd2 is None else jd2
        for jd in np.arange(jd1, jd2):
            print " --- %s --- " % pl.num2date(jd).strftime('%Y-%m-%d')
            filename = os.path.join(
                self.datadir, self.generate_filename(jd,fld,fldtype) + ".nc")
            if delall:
                for fn in glob.glob(filename + "*"):
                    print "Deleted %s" % fn
                    os.remove(fn)
            print "Checking files"
            if not os.path.isfile(filename[:-3] + '.npz'):
                try:
                    self.load(fld, fldtype, jd=jd, verbose=True)
                except IOError:
                    print "Downloading failed. Trying to remove old files."
                    try:
                        os.remove(filename)
                    except:
                        pass
                    try:
                        self.load(fld,fldtype,jd=jd,verbose=True)
                    except:
                        print ("   ###   Warning! Failed to add %s   ###" %
                               os.path.basename(filename))
                print "\n"
            else:
                print "found"

             
    def load(self, fld, fldtype="DAY", nan="nan", **kwargs):
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
        self.filename = self.generate_filename(fld, fldtype)
        self.vprint( "datadir:        %s" % self.datadir)
        self.vprint( "filename:       %s" % os.path.basename(self.filename))
        self._l3read(fld,nan=nan)

    @property
    def landmask(self):
        """Return a landmask"""
        if not hasattr(self,'_landmask'):
            self.load('par','CU', nan="nan")
            self._landmask = np.isnan(self.par)
        return self._landmask
    
    def add_mnclim(self):
        """Add a list with file name dates for Monthly climatologies"""
        datelist = self._retrieve_datestamps(self.a_mc_url_9km)
        self.mc_datedict = {}
        for dstr in datelist[1:]:
            mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
            self.mc_datedict[mn] = dstr

    def add_filepreflist(self, fldtype="MC"):
        """Add a list with file name dates for Monthly climatologies"""
        if hasattr(self, "%s_fileprefs" % fldtype):
            return

        predict = {}
        if "mc" in fldtype.lower():
            datelist = self._retrieve_datestamps(self.a_mc_url_9km)
            self.mc_datedict = {}
            for dstr in datelist[1:]:
                mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
                self.mc_datedict[mn] = dstr
        elif "mo" in fldtype.lower():
            datelist = self._retrieve_datestamps(self.a_mo_url_9km)
            for dstr in datelist[1:]:
                mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
                yr = int(dstr[:4])
                predict[yr*100 + mn] = dstr
        setattr(self, "%s_fileprefs" % fldtype, predict)

    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['KD490_Kd_490',         'km',  0,      100],
                   'chl': ['CHL_chlor_a',          'km',  0.001,  200],
                   'poc': ['POC_poc',              'km',  0,     1000],
                   'bbp': ['GIOP01_bbp_443_giop',  'km',  0,    10000],
                   'sst': ['SST_sst',              'km', -2,       45],
                   'par': ['PAR_par',              'km', -5,      200],
                   'flh': ['FLH_nflh',             'km',  0,    10000],
                   'ipar': ['FLH_ipar',            'km',  0,    10000],
                   'eup': ['KDLEE_Zeu_lee',        'km',  0,    10000],
                   'eup_lee': ['KDLEE_Zeu_lee',    'km',  0,    10000],
                   'eup_mor': ['KDMOREL_Zeu_morel','km',  0,    10000]
                   }

    def download(self, filename):
        """Download a missing file from GSFC's website"""
        uri = "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/"
        url = "%s%s" % (uri, os.path.basename(filename))
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.retrive_file(url, local_filename=filename)
        """
        try:
            response = urllib2.urlopen(url)
        except:
            raise IOError, "File not found on the server.\n tried %s" % url
        output = open(filename, 'wb')
        output.write(response.read())
        output.close()
        """
    def _l3read(self,fld='',nan="nan"):
        """ Read a L3 mapped file and add field to current instance"""
        self.filename = self.datadir + '/' + self.filename
        if not fld: fld = self.filename.rsplit('_')[-2]

        if os.path.isfile(self.filename[:-3] + '.npz') & (self.use_npz == True):
            self.vprint( "Found npz file.")
            l3m_data,base,intercept,slope = self._l3read_npz(fld)
        else:
            self.vprint( "No npz file.")
            if not os.path.isfile(self.filename):
                if getattr(self, "no_download", False):
                    raise IOError, "File Missing"
                else:
                    print "File missing, downloading from GSFC"
                    self.download(self.filename)

            l3m_data = self._l3read_nc4()
            base=-999; intercept=0; slope=1

        if base != -999:
            field = self.__dict__[fld] = base**((slope * l3m_data) + intercept)
        else:
            field = self.__dict__[fld] = ((slope * l3m_data) + intercept)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mask = (field < self.minval) | (field > self.maxval)
            if nan=="nan":
                field[mask] = np.nan

        if self.use_npz == True:
            if os.path.isfile(self.filename[:-3] + '.npz'):
                self.vprint( "npz file already exists, no need to save.")
            elif (self.llat.shape) != (self.jmt,self.imt):
                self.vprint( "Not global grid, can't generate npz file.")
            else:
                self.add_ij()
                self._l3write_npz(l3m_data[~mask], self.imat[~mask],
                                  self.jmat[~mask], base, intercept, slope)
            if (self.only_npz == True) & os.path.isfile(self.filename):
                os.remove(self.filename)
        else:
            self.vprint( "use_npz set to False.")
            
            
    def _l3read_nc4(self):
        self.vprint( "Reading netCDF4 file")
        print(self.filename)
        nc = Dataset(self.filename)
        nc.set_auto_mask(False)
        nc.set_auto_scale(True)
        
        var         = nc.variables[nc.variables.keys()[0]]
        field       = var[self.j1:self.j2, self.i1:self.i2].copy()
        try:
            self.minval = var.valid_min
        except AttributeError:
            self.minval = var.display_min
        try:
            valid_max = var.valid_max
        except AttributeError:
            valid_max = 0
        try:
            display_max = var.display_max
        except AttributeError:
            display_max = 0
        self.maxval = max(valid_max, display_max)
            
        start_jd    = pl.datestr2num(nc.time_coverage_start)
        end_jd      = pl.datestr2num(nc.time_coverage_end)
        self.jd     = ((start_jd + end_jd)/2)
        self.date   = pl.num2date(self.jd)
        return field

    def _l3read_npz(self, fld):
        """Read an npz file with field data as vectors"""
        self.vprint( "Reading npz file")
        l3m_data = self.llat * 0 + self.vc[fld][2] - 1
        fH = np.load(self.filename[:-3] + '.npz')
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
        self.vprint( "Writing npz file")
        kwargs = {'dvec' : dvec,
                  'ivec' : ivec.astype(np.int16),
                  'jvec' : jvec.astype(np.int16),
                  'base' : base, 'intercept':intercept, 'slope':slope,
                  'jd'   : self.jd, 'date':self.date}
        np.savez_compressed(self.filename[:-3] + '.npz', **kwargs)



class MODIS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "A"
        super(MODIS, self).__init__(**kwargs)

class SeaWIFS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "S"
        super(SeaWIFS, self).__init__(**kwargs)

class OCTS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "O"
        super(OCTS, self).__init__(**kwargs)

class CZCS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "C"
        super(CZCS, self).__init__(**kwargs)

    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['MO_KD490_Kd_490',         'km',  0,      100],
                   'chl': ['CHL_chlor_a',       'km',  0.001,  200],
                   'poc': ['MO_POC_poc',           'km',  0,     1000],
                   'bbp': ['GIOP01_bbp_443_giop',  'km',  0,    10000],
                   'sst': ['MO_SST_sst',           'km', -2,       45],
                   'par': ['MO_PAR_par',           'km', -5,      200],
                   'flh': ['FLH_nflh',             'km',  0,    10000],
                   'ipar': ['FLH_ipar',            'km',  0,    10000],
                   'eup': ['KDLEE_Zeu_lee',        'km',  0,    10000],
                   'eup_lee': ['KDLEE_Zeu_lee',    'km',  0,    10000],
                   'eup_mor': ['KDMOREL_Zeu_morel','km',  0,    10000]
                   }

class VIIRS(Base):
    def __init__(self, res="9km", **kwargs):
        self.res = res        
        self.fp = "V"
        super(VIIRS, self).__init__(**kwargs)

    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['NPP_KD490_Kd_490',     'km',  0,      100],
                   'chl': ['NPP_CHL_chlor_a',      'km',  0.001,  200],
                   'poc': ['NPP_POC_poc',          'km',  0,     1000],
                   'bbp': ['GIOP01_bbp_443_giop',  'km',  0,    10000],
                   'sst': ['SST_sst',              'km', -2,       45],
                   'par': ['PAR_par',              'km', -5,      200],
                   'flh': ['FLH_nflh',             'km',  0,    10000],
                   'ipar': ['FLH_ipar',            'km',  0,    10000],
                   'eup': ['KDLEE_Zeu_lee',        'km',  0,    10000],
                   'eup_lee': ['KDLEE_Zeu_lee',    'km',  0,    10000],
                   'eup_mor': ['KDMOREL_Zeu_morel','km',  0,    10000]
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
            print jd,pl.num2date(jd).year,mn
        return MCmin
