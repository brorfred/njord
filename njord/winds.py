import os, os.path
from datetime import datetime as dtm
import ftplib
import urlparse

import scipy
import numpy as np
import pylab as pl
import netCDF4

import bs4
import requests

import base
import gmtgrid

class Seawinds(base.Grid):
    """Read jpl Seawinds fields"""
    def __init__(self, **kwargs):
        """Initialize the class with stuff from base.Grid"""
	super(Seawinds, self).__init__(**kwargs)
	
    def setup_grid(self):
	"""Setup lat-lon matrices for Seawinds"""
        if not os.path.isfile(self.gridfile):
            self.retrive_file(self.dataurl, self.gridfile)
        try:
            n = netCDF4.Dataset(self.gridfile, 'r')
        except:
            print 'Error opening the gridfile %s' % self.gridfile
            raise
        self.lat = n.variables['lat'][:]
        self.gmt = gmtgrid.Shift(n.variables['lon'][:].copy())
        self.lon = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)

    def load(self, fld="nwnd", **kwargs):
	"""Load field for a given julian date. Returns u,v, or nwnd(windnorm)"""
        self._timeparams(**kwargs)
        filename = os.path.join(self.datadir,
                            "uv%04i%02i%02i.nc" % (self.yr, self.mn, self.dy))
        if not os.path.isfile(filename):
            self.download(filename)
        try:
            nc = netCDF4.Dataset(filename)
        except:
            os.remove(filename)
            self.download(filename)
            try:
                nc = netCDF4.Dataset(filename)
            except:
                filename = filename.rstrip(".nc") + "rt.nc"
                if not os.path.isfile(filename):
                    self.download(filename)
                try:
                    nc = netCDF4.Dataset(filename)
                except TypeError:
                    os.remove(filename)
                    self.download(filename)
                    nc = netCDF4.Dataset(filename)
                    
        u = nc.variables['u'][:].copy()
        v = nc.variables['v'][:].copy()
        u[u<-999] = np.nan
        v[v<-999] = np.nan
        if (fld=="u") | (fld=="uvel"):
            self.uvel = self.gmt.field(np.squeeze(u))
        elif (fld=="v") | (fld=="vvel"):
            self.vvel = self.gmt.field(np.squeeze(v))
        else:
            self.nwnd = self.gmt.field(np.squeeze(np.sqrt(u**2 + v**2)))

    def download(self, filename):
        try:
            self.retrive_file(self.dataurl, filename)
        except ftplib.error_perm:
            return False


class CCMP(base.Grid):
    """Read jpl CCMP fields
    http://podaac.jpl.nasa.gov/dataset/CCMP_MEASURES_ATLAS_L4_OW_L3_0_WIND_VECTORS_FLK
     ftp://podaac-ftp.jpl.nasa.gov/allData/ccmp/

    """
    def __init__(self, **kwargs):
        """Initialize the class with stuff from base.Grid"""
	super(CCMP, self).__init__(**kwargs)
	
    def setup_grid(self):
	"""Setup lat-lon matrices for CCMP"""
	try:
	    gc = netCDF4.Dataset(self.gridfile, 'r')
        except:
            print 'Error opening the gridfile %s' % datadir + filename
            raise
        self.lat = gc.variables['lat'][:]
	self.gmt = gmtgrid.Shift(gc.variables['lon'][:].copy())
        self.lon = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)

    def load(self, fld="nwnd", **kwargs):
	"""Load field for a given julian date. Returns u,v, or nwnd(windnorm)"""
        self._timeparams(**kwargs)
	filename = os.path.join(self.datadir,
				"analysis_%04i%02i%02i_v11l30flk.nc" %
                                  	(self.yr,self.mn,self.dy))
        if os.path.isfile(filename):
            nc = netCDF4.Dataset(filename)
        else:
            raise IOError, 'Error opening the windfile %s' % filename
	uH = nc.variables['uwnd']
	vH = nc.variables['vwnd']
	uvel = self.gmt.field(uH.data.copy()) * uH.scale_factor
	vvel = self.gmt.field(vH.data.copy()) * vH.scale_factor
	
	uvel[uvel<(uH.missing_value * uH.scale_factor)] = np.nan
	vvel[vvel<(vH.missing_value * vH.scale_factor)] = np.nan
	if (fld=="u") | (fld=="uvel"):
	    self.uvel = gmtgrid.convert(np.squeeze(u), self.gr)
	elif (fld=="v") | (fld=="vvel"):
	    self.vvel = gmtgrid.convert(np.squeeze(v), self.gr)
	else:
   	    self.nwnd = gmtgrid.convert(np.squeeze(np.sqrt(u**2 + v**2)),self.gr)




class ncep:

    def __init__(self,datadir = "/projData/ncep/"):
        filename = "vwnd.sig995.2008.nc"        
        try:
            n = pycdf.CDF(datadir + filename)
        except:
            print 'Error opening the gridfile %s' % datadir + filename
            raise
        self.lat = n.var('lat')[:]
        self.lon,self.gr = gmtgrid.config(n.var('lon')[:],0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.datadir = datadir

    def load(self,jd):
        yr = pl.num2date(jd).year
        yd = int((jd - pl.date2num(dtm(yr,1,1))) * 4)
        ufile = "uwnd.sig995.%04i.nc" % yr
        vfile = "vwnd.sig995.%04i.nc" % yr
        try:
            un = pycdf.CDF(self.datadir + ufile)
        except:
            print 'Error opening the windfile %s' % datadir + ufile
            raise
        try:
            vn = pycdf.CDF(self.datadir + vfile)
        except:
            print 'Error opening the windfile %s' % datadir + vfile
            raise    
        u = un.var('uwnd')[yd,:,:] * 0.01 + 225.45
        v = vn.var('vwnd')[yd,:,:] * 0.01 + 225.45
        nwnd = gmtgrid.convert(np.sqrt(u**2 + v**2),self.gr)
        nwnd[nwnd>200]=np.nan
        return nwnd

class Quikscat(base.Grid):
    
    def __init__(self, **kwargs):
        super(Quikscat, self).__init__(**kwargs)

    def setup_grid(self):
        self.latvec = np.linspace(-89.875, 89.875, 720)
        self.gmt = gmtgrid.Shift(np.linspace(0, 360, 1440))
        self.lonvec = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)
        
    def load(self, fldname, **kwargs):

        def merge(fldname, fldname1, fldname2):
            fld1 = self.get_field(fldname1, **kwargs)
            fld2 = self.get_field(fldname2, **kwargs)
            setattr(self, fldname, np.nanmedian((fld1, fld2),axis=0))
        
        if fldname == "uvel":
            merge("uvel", "des_avg_wind_vel_u", "asc_avg_wind_vel_u")
            return
        if fldname == "vvel":
            merge("vvel", "des_avg_wind_vel_v", "asc_avg_wind_vel_v")
            return
        if fldname == "nwnd":
            uvel = self.get_field("uvel", **kwargs)
            vvel = self.get_field("vvel", **kwargs)
            setattr(self, "nwnd", np.sqrt(uvel**2 + vvel**2))
            return
        
        self._timeparams(**kwargs)
        filename = self.filedict(self.yr)[self.jd] + ".gz.nc"
        fn = os.path.join(self.datadir,  filename)
        if not os.path.isfile(fn):
            url = "%s/%s/%s" % (self.dataurl, self.yr, filename)
            self.retrive_file(url, fn)
        fld = self.gmt.field(netCDF4.Dataset(fn).variables[fldname][:].copy())
        fld = fld[self.j1:self.j2, self.i1:self.i2]
        fld[fld==0] = np.nan
        setattr(self, fldname, fld)

    def filedict(self, year):
        filename = os.path.join(self.datadir,  "fdict_%04i.npz" % year)
        if os.path.isfile(filename):
            return np.load(filename)["fdict"].item()
        else:
            req = requests.get("%s/%04i/%s" %(self.dataurl,year,self.flistpage))
            soup = bs4.BeautifulSoup(req.text,"html.parser")
            taglist = soup.find_all("a",itemprop="contentUrl", string="html")
            fdict = {}
            for tag in taglist:
                stamp = ".".join(tag.attrs['href'].split(".")[:2])
                datestr = stamp.split("_")[-1].split(".")[0]
                jd = int(pl.datestr2num("%s-01-01" % datestr[:4]) +
                         int(datestr[4:]))
                fdict[jd] = stamp
            np.savez(filename, fdict=fdict)
            return fdict
            
class CORE2:

    def __init__(self,datadir = "/projData/CORE2/"):
        self.jdbase = pl.date2num(dtm(1948,1,1))+15
        self.datadir = datadir
        filename = "u_10.2005.05APR2010.nc"
        try:
            n = netCDF4.Dataset(datadir + filename)
        except:
            print 'Error opening the gridfile %s' % datadir + filename
            raise
        self.lat = n.variables['LAT'][:]
        self.gmt = gmtgrid.Shift(n.variables['LON'][:].copy())
        self.lon = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        n.close()
	
    def load(self,jd1,jd2=None):
        yr = pl.num2date(jd1).year
        mn = pl.num2date(jd1).month
        dy = pl.num2date(jd1).day
        filesuff = ".%04i.05APR2010.nc" %(yr)
        try:
            nu = netCDF4.Dataset(self.datadir + "u_10" + filesuff)
            nv = netCDF4.Dataset(self.datadir + "v_10" + filesuff)
        except:
            print 'Error opening the windfile %s%s' % (self.datadir,filesuff)
            raise
        jdvec = nu.variables['TIME'][:] + self.jdbase
        t1 = int(np.nonzero(jdvec>jd1)[0].min())
        if not jd2:
            jd2 = jd1+1
            t2 = t1+4
        elif jd2<jdvec.max():
            t2 = int(t1+4*(jd2-jd1))
        else:
            jd11 = jd1
            jd22 = jdvec.max() -1
            while np.floor(jd2)>=np.ceil(jd22):
                print jd22-jd11,jd2-jd11,jd11,jd22
                wnd1 = self.load(jd11,jd22)
                jd11 = jd22 + 1
                jd22 += 365

            return
        djd = np.ceil(jd2 - jd1)
        wndjd = np.zeros( ( (djd,) + self.llon.shape) ) 
        self.uwnd = self.gmt.field(nu.variables['U_10_MOD'][t1:t2,...].copy())
        self.vwnd = self.gmt.field(nv.variables['V_10_MOD'][t1:t2,...].copy())
        nwnd = np.sqrt(self.uwnd**2 + self.vwnd**2)
        for t in np.arange(0,djd*4,4):
            wndjd[t/4,...] = np.mean(nwnd[t:t+4,...],axis=0)
        nu.close()
        nv.close()
        self.nwnd = np.squeeze(wndjd)
