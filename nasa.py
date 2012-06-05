import os, os.path
from datetime import datetime as dtm
import urllib
import urllib2

import numpy as np
import pylab as pl
import numpy ,pyhdf ,os.path
from scipy.spatial import KDTree, cKDTree

from pyhdf.SD import SD, SDC
import njord

class MODIS(njord.Njord):

    def __init__(self ,res="9km", ijarea=None,i1=None,i2=None,j1=None,j2=None,
                 latlon=None,lat1=None,lat2=None,lon1=None,lon2=None):

        if res == "9km":
            fi1,fi2,fj1,fj2 = (0000 ,2160, 0, 4320)
            incr = 360.0/4320.0
            self.datadir = "/projData/sat/A9km/"
        elif res == "4km":
            fi1,fi2,fj1,fj2 = (0000 ,4320, 0, 8640)
            incr = 360.0/8640.0
            self.datadir = "/projData/sat/A4km/"
        iR    = numpy.arange(fj1, fj2)
        jR    = numpy.arange(fi1, fi2)
        [x,y] = numpy.meshgrid(iR,jR)
        self.llat = (  90 - y*incr - incr/2)
        self.llon = (-180 + x*incr + incr/2)
        i1,i2,j1,j2 = (fi1, fi2, fj1, fj2)
     
        if ijarea is not None:
            i1,i2,j1,j2 = ijarea
        elif latlon is not None:
            i1,i2,j1,j2 = self.llrect2ij(*latlon)   
        if lat1 is not None:
            i2 = int(np.nonzero(self.llat[:,0] < lat1)[0].min())
        if lat2 is not None:
            i1 = int(np.nonzero(self.llat[:,0] > lat2)[0].max())
        if lon1 is not None:
            j2 = int(np.nonzero(self.llon[0,:] > lon2)[0].min())
        if lon2 is not None:
            j1 = int(np.nonzero(self.llon[0,:] < lon1)[0].max())
        self.i1,self.i2,self.j1,self.j2 = i1,i2,j1,j2

        self.llat = (  90 - y*incr - incr/2)[i1:i2,j1:j2]
        self.llon = (-180 + x*incr + incr/2)[i1:i2,j1:j2]
        self.lat = self.llat[:,0]
        self.lon = self.llon[0,:]
        self.res = res
        self.add_vc()

    def add_landmask(self):
        """Add a landmask field to the current instance"""
        if self.res == "4km":
            file = 'A20021852011273.L3m_CU_PAR_par_4km'
            h = SD(self.datadir + file, SDC.READ)
            fld = h.select('l3m_data')[self.i1:self.i2, self.j1:self.j2]
            fld[fld!=-32767.] = 0
            fld[fld==-32767.] = 1
            self.landmask = fld.astype(np.int8)
        elif self.res == "9km":
            file = 'A20021852011031.L3m_CU_CHL_chlor_a_9km'
            h = SD(self.datadir + file, SDC.READ)
            fld = h.select('l3m_data')[self.i1:self.i2, self.j1:self.j2]
            fld[fld!=-32767.] = 0
            fld[fld==-32767.] = 1
            self.landmask = fld.astype(np.int8)
        else:
            print "Unknown Resolution"
            raise

    def add_mnclim(self):
        """Add a list with file name dates for Monthly climatologies"""
        self.mnclim = [
            "A20030012010031",
            "A20030322010059",
            "A20030602010090",
            "A20030912010120",
            "A20031212010151",
            "A20031522010181",
            "A20021822010212",
            "A20022132010243",
            "A20022442010273",
            "A20022742010304",
            "A20023052009334",
            "A20023352009365"]
        
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

    def l2read(self ,filename,fld):
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

    def l3read(self ,filename,fld=''):
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
            sd.select('l3m_data')[self.i1:self.i2, self.j1:self.j2],65535)
        Intercept = sd.attributes()['Intercept']
        Slope     = sd.attributes()['Slope']
        try:    
            Base  = sd.attributes()['Base']
            self.__dict__[fld] = Base**((Slope*l3m_data) + Intercept)
        except KeyError:
            self.__dict__[fld] = ((Slope*l3m_data) + Intercept)
        self.__dict__[fld][self.__dict__[fld] < self.vc[fld][2]] = np.nan

        tmp = self.__dict__[fld]
        tmp = tmp[~np.isnan(tmp)]
        if len(tmp)>0:
            print  tmp.min(), tmp.max()
            
        start_iso = (pl.date2num(dtm(
                    sd.attributes()['Period Start Year'],1,1)) + 
                     sd.attributes()['Period Start Day'] - 1)
        end_iso   = (pl.date2num(dtm(
                    sd.attributes()['Period End Year'],1,1)) + 
                     sd.attributes()['Period End Day'] - 1)
        self.date = pl.num2date((start_iso+end_iso)/2)
        self.iso  = ((start_iso+end_iso)/2)

        if zipped:
            try:
                os.system('pbzip2 ' + filename)
            except:
                raise IOerror( "Compression of " + filename + " failed.")

    def load(self,fld,fldtype="DAY",yr=2005,yd=300,mn=1,dy=0,jd=0):
        """Load a satellite field."""
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
            filename = self.mnclim[mn-1]
        elif fldtype == "DAY":
            datestr = "%i%03i" % (yr, yd)
        elif fldtype == "8D":
            yd1 = np.arange(1,365,8)
            yd2 = np.arange(8,370,8)
            yd2[-1] = ydmax
            pos = np.nonzero(yd >= yd1)[0].max()
            datestr = "%i%03i%i%03i" % (yr, yd1[pos], yr, yd2[pos])
        elif fldtype == "CU":
            datestr = "A20021852011059"
        else:
            print "Field type not included"
        filename = ("A%s.L3m_%s_%s_%s%s" % (datestr, fldtype,
                                            self.vc[fld][0],
                                            self.res[0],
                                            self.vc[fld][1]))
        self.l3read(filename,fld)



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
