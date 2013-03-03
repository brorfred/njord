import os.path
from datetime import datetime as dtm

import scipy
import numpy as np
import pylab as pl
from scipy.stats import nanmean
from scipy.io import netcdf_file

import base
import gmtgrid
import reuerflux
import bln


class Winds(object):
    """ Meta class for winds """
    def timeseries(self, jd1,jd2,i,j):
        """Generate a timeserie for a gridcell"""
        vec = np.array([self.load(jd)[i,j] for jd in np.arange(jd1,jd2+1)])
        return vec    

class Seawinds(base.Grid):
    """Read jpl Seawinds fields"""
    def __init__(self, **kwargs):
        """Initialize the class with stuff from base.Grid"""
	super(Seawinds, self).__init__(**kwargs)
	
    def setup_grid(self):
	"""Setup lat-lon matrices for Seawinds"""
        try:
	    n = netcdf_file(self.gridfile, 'r')
        except:
            print 'Error opening the gridfile %s' % datadir + filename
            raise
        self.lat = n.variables['lat'][:]
        self.lon,self.gr = gmtgrid.config(n.variables['lon'][:].copy(), 0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)

    def load(self, fld="nwnd", jd=733805.0):
	"""Load field for a given julian date. Returns u,v, or nwnd(windnorm)"""
        yr = pl.num2date(jd).year
        mn = pl.num2date(jd).month
        dy = pl.num2date(jd).day
        filename1 = "uv%04i%02i%02i.nc" %(yr,mn,dy)
        filename2 = "uv%04i%02i%02irt.nc" %(yr,mn,dy)
        if os.path.isfile(datadir + filename1):
            nc = netcdf_file(datadir + filename1)
        elif os.path.isfile(datadir + filename2):
            nc = netcdf_file(datadir + filename2)
        else:
            raise IOError, 'Error opening the windfile %s' % datadir + filename1
            raise IOError, 'Error opening the windfile %s' % datadir + filename2
	u = nc.variables['u'][:].copy()
        v = nc.variables['v'][:].copy()
	u[u<-999] = np.nan
	v[v<-999] = np.nan
	if (fld=="u") | (fld=="uvel"):
	    self.uvel = gmtgrid.convert(np.squeeze(u), self.gr)
	elif (fld=="v") | (fld=="vvel"):
	    self.vvel = gmtgrid.convert(np.squeeze(v), self.gr)
	else:
   	    self.nwnd = gmtgrid.convert(np.squeeze(np.sqrt(u**2 + v**2)),self.gr)



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
	    gc = netcdf_file(self.gridfile, 'r')
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
            nc = netcdf_file(filename)
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

class Quikscat(Winds):
    
    def __init__(self,datadir = "/projData/QUIKSCAT/"):
        self.lat,self.lon = bln.grid()
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.datadir = datadir
        
    def load(self,jd):
        u,v = bln.readuv_day(jd)
        nwnd =np.sqrt(u**2 + v**2)
        #nwnd[nwnd>200]=np.nan
        return nwnd

class CORE2:

    def __init__(self,datadir = "/projData/CORE2/"):
        self.jdbase = pl.date2num(dtm(1948,1,1))+15
        self.datadir = datadir
        try:
            n = pycdf.CDF(datadir + "u_10.2005.05APR2010.nc")
        except:
            print 'Error opening the gridfile %s' % datadir + filename
            raise
        self.lat = n.var('LAT')[:]
        self.lon,self.gr = gmtgrid.config(n.var('LON')[:],dim=0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        n.close()
	
    def load(self,jd1,jd2=None):
        yr = pl.num2date(jd1).year
        mn = pl.num2date(jd1).month
        dy = pl.num2date(jd1).day
        filesuff = ".%04i.05APR2010.nc" %(yr)
        try:
            nu = pycdf.CDF(self.datadir + "u_10" + filesuff)
            nv = pycdf.CDF(self.datadir + "v_10" + filesuff)
        except:
            print 'Error opening the windfile %s' % datadir + "*" + filesuff
            raise
        jdvec = nu.var('TIME')[:] + self.jdbase
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
        u = nu.var('U_10_MOD')[t1:t2,...]
        v = nv.var('V_10_MOD')[t1:t2,...]
        nwnd = gmtgrid.convert(np.sqrt(u**2 + v**2), self.gr)
        for t in np.arange(0,djd*4,4):
            wndjd[t/4,...] = np.mean(nwnd[t:t+4,...],axis=0)

        #nwnd[nwnd>200]=np.nan
        nu.close()
        nv.close()
        return np.squeeze(wndjd)
