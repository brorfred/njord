import os, os.path
import glob
import subprocess as sbp
from datetime import datetime as dtm
import urllib
import urllib2
import re
from distutils import spawn

import numpy as np
import pylab as pl

from pyhdf.SD import SD, SDC
import base
from utils import yrday

if spawn.find_executable("pbzip2"):
    ZIPCMD = "pbzip2"
elif spawn.find_executable("bzip2"):
    ZIPCMD = "bzip2"
else:
    raise IOError, "Couldn't find a bzip2 executable. Needed to unzip files"

class Base(base.Grid):

    def __init__(self, **kwargs):
        super(Base, self).__init__(**kwargs)
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.add_vc()
        self.datadir = '%s/%s%s/' % (self.datadir, self.fp, self.res) 
        self.only_npz = kwargs.get('only_npz', False)

    def setup_grid(self):
        """Create matrices with latitudes and longitudes for the t-coords"""
        if self.res == "9km":
            i1,i2,j1,j2 = (0000 ,4320, 0, 2160)
        elif self.res == "4km":
            i1,i2,j1,j2 = (0000 ,8640, 0, 4320)
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

    def generate_filename(self, jd, fld='chl', fldtype="DAY"):
        """Generate filename"""
        self._timeparams(jd=jd)
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
        return("%s%s.L3m_%s_%s_%s%s" % (self.fp, datestr, fldtype,
										self.vc[fld][0], self.res[0],
                                        self.vc[fld][1]))

    def refresh(self, fld="chl", fldtype="DAY", jd1=None, jd2=None, delall=False):
        """ Read a L3 mapped file and add field to current instance"""
        jd1 = pl.datestr2num('2003-01-01') if jd1 is None else jd1
        jd2 = int(pl.date2num(dtm.now())) - 1  if jd2 is None else jd2
        for jd in np.arange(jd1, jd2):
            print " --- %s --- " % pl.num2date(jd).strftime('%Y-%m-%d')
            filename = os.path.join(self.datadir,
                                    self.generate_filename(jd,fld,fldtype))
            if delall:
                for fn in glob.glob(filename + "*"):
                    print "Deleted %s" % fn
                    os.remove(fn)
            print "Checking %s" % filename + '.npz'
            if not os.path.isfile(filename + '.npz'):
                try:
                    self.load(fld, fldtype, jd=jd, verbose=True)
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
        self.filename = self.generate_filename(self.jd, fld, fldtype)
        self.vprint( "Filename is %s" % (self.filename))
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
        url = "http://%s/%s" % (self.gsfc_host, self.a_mc_url_9km)
        datelist = []
        for line in urllib.urlopen(url):
            if "cgi/getfile" in line:
                datelist.append(re.findall(r'getfile/%s([^E]+).L3m' % self.fp,
                                           line)[0][:14])
        self.mc_datedict = {}
        for dstr in datelist[1:]:
            mn,_ = yrday.mndy(int(dstr[:4]),int(dstr[4:7]))
            self.mc_datedict[mn] = dstr
        
    def add_vc(self):
        """Add a dict with filename variable components"""
        self.vc = {'k49': ['KD490_Kd_490',         'km',  0,      100],
                   'chl': ['CHL_chlor_a',          'km',  0.001,  200],
                   'poc': ['POC_poc',              'km',  0,     1000],
                   'bbp': ['GIOP01_bbp_443_giop',  'km',  0,    10000],
                   'sst': ['SST',                  '',   -2,       45],
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
        url = "%s%s.bz2" % (uri, os.path.basename(filename))
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
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
            self.vprint( "Found npz file.")
            l3m_data,base,intercept,slope = self._l3read_npz(fld)
            zipped = False
        else:
            self.vprint( "Didn't find a npz file.")
            zipfname = self.filename + '.bz2'
            if not (os.path.isfile(self.filename) | 
                    os.path.isfile(zipfname)):
                print "File missing, downloading from GSFC"
                self.download(self.filename)
            self.vprint( "Found bz2 file: %s" % zipfname)
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

        if os.path.isfile(self.filename + '.npz'):
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
            self._try_to_zip(zipped)


    def _l3read_hdf(self, fieldname='l3m_data'):
        self.vprint( "Reading hdf file")
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
        self.vprint( "Reading npz file")
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
        self.vprint( "Writing npz file")
        kwargs = {'dvec' : dvec,
                  'ivec' : ivec.astype(np.int16),
                  'jvec' : jvec.astype(np.int16),
                  'base' : base, 'intercept':intercept, 'slope':slope,
                  'jd'   : self.jd, 'date':self.date}
        np.savez_compressed(self.filename + '.npz', **kwargs)

    def _try_to_unzip(self, zipfname):
        """Unzip file if exists and and is a valid bzip2 file"""
        if not os.path.isfile(self.filename):
            self.vprint( "Uncompressing file")
            try:
                err = sbp.call([ZIPCMD, "-d", zipfname])
            except OSError as det:
                raise OSError, "Error when unzipping %s: %s" % (zipfname,det)
                
            if err == 1:
                print "Decompression of " + zipfname + " failed."
                print "Trying to download again"
                self.download(self.filename)
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
