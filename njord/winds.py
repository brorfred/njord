import os, os.path
from datetime import datetime as dtm
import ftplib
import gzip

import scipy
import numpy as np
import njord.utils.mpl_dates as pl
import netCDF4

import bs4
import requests

from . import base
from njord.utils import gmtgrid

class Seawinds(base.Grid):
    """Read jpl Seawinds fields"""
    def __init__(self, **kwargs):
        """Initialize the class with stuff from base.Grid"""
        super(Seawinds, self).__init__(**kwargs)
        self.fieldlist = ["uwnd", "uwnd", "nwnd"]
        if not hasattr(self, "maxjd"):
            self.maxjd = int(pl.date2num(dtm.now())) - 3

    def setup_grid(self):
        """Setup lat-lon matrices for Seawinds"""
        if not os.path.isfile(self.gridfile):
            self.retrive_file(self.dataurl, self.gridfile)
        try:
            n = netCDF4.Dataset(self.gridfile, 'r')
        except:
            raise IOError("Cant't open the gridfile %s" % self.gridfile)
        self.latvec = n.variables['lat'][:]
        self.gmt = gmtgrid.Shift(n.variables['lon'][:].copy())
        self.lonvec = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    def load(self, fld="nwnd", **kwargs):
        """Load wind field"""
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

        def get_fld(varobj):
            raw = varobj[:]
            fld = raw.data
            fld[raw.mask] = np.nan
            return self.gmt.field(np.squeeze(fld))[self.ijslice]
        
        if (fld=="u") | (fld=="uvel") | (fld=="uwnd"):
            setattr(self, "uwnd", get_fld(nc.variables['u']))
        elif (fld=="v") | (fld=="vvel") | (fld=="uwnd"):
            setattr(self, "vwnd", get_fld(nc.variables['v']))
        else:
            uwnd = get_fld(nc.variables['u'])            
            vwnd = get_fld(nc.variables['v'])        
            setattr(self, "nwnd", np.sqrt(uwnd**2 + vwnd**2))

    def download(self, filename):
        """Download missing data file from server"""
        try:
            self.retrive_file(self.dataurl, filename)
        except ftplib.error_perm:
            return False

    @property
    def landmask(self):
        if not hasattr(self,"_landmask"):
            filename = os.path.join(self.datadir, "topo15g.asc")
            if not os.path.isfile(filename):
                url = "/".join(self.dataurl.split("/")[:-5])+"/"
                self.retrive_file(url, filename)
            self._landmask = np.loadtxt(filename) > 0
            if len(self._landmask) == 1:
                raise IOError("%s empty, try to delete the file and retry.")
        return self.gmt.field(np.squeeze(self._landmask))[self.ijslice]

    
class CCMP(base.Grid):
    """Read jpl CCMP fields"""
    def __init__(self, **kwargs):
        """Initialize the class with stuff from base.Grid"""
        super(CCMP, self).__init__(**kwargs)
        self.fieldlist = ["uwnd", "uwnd", "nwnd"]

    def setup_grid(self):
        """Setup lat-lon matrices for CCMP"""
        if not os.path.isfile(self.gridfile):
            self.load("uwnd", yr=2007, mn=12, dy=1)
        gc = netCDF4.Dataset(self.gridfile, 'r')
        self.latvec = gc.variables['lat'][:]
        self.gmt = gmtgrid.Shift(gc.variables['lon'][:].copy())
        self.lonvec = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    def download(self, filename):
        url = "%s/%04i/%02i/" % (self.dataurl, self.yr, self.mn)
        self.retrive_file(url, self.filename + ".gz")
        with gzip.open(self.filename + ".gz", 'rb') as gzfile:
                with open(self.filename, 'wb') as ncfile:
                    ncfile.write( gzfile.read() )

    def load(self, fldname="nwnd", use_mask=False, **kwargs):
        """Load wind field"""
        if fldname == "nwnd":
            uwnd = self.get_field("uwnd", **kwargs)
            vwnd = self.get_field("uwnd", **kwargs)
            setattr(self, "nwnd", np.sqrt(uwnd**2 + vwnd**2))
            return
        
        self._timeparams(**kwargs)
        filename = ("analysis_%04i%02i%02i_v11l30flk.nc" %
            (self.yr,self.mn,self.dy))
        self.filename = os.path.join(self.datadir, filename)
        if not os.path.isfile(filename):
            self.download(self.filename)
        raw = self.gmt.field(
        netCDF4.Dataset(self.filename).variables[fldname][:])
        fld = raw.data
        if use_mask:
            fld[raw.mask] = np.nan
        if "hr" in kwargs:
            fld = fld[kwargs["hr"]/6,...][self.ijslice]
        else:
            fld = np.nanmean(fld, axis=0)[self.ijslice]
        setattr(self, fldname, fld)
        
    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            if not os.path.isfile(self.landfile):
                self.retrive_file(self.landurl, self.landfile)
            mask = np.load(self.landfile)["aggr_nobs"]
            self._landmask = ~mask.astype(bool)[self.ijslice]
        return self._landmask
            
class NCEP(base.Grid):

    def __init__(self, version=1, **kwargs):
        super(NCEP, self).__init__(**kwargs)
        if version not in [1,]:
            raise TypeError("Only version 1 works for the moment")
        setattr(self, "datadir", os.path.join(self.datadir, "v%i" % version))
        if version == 2:
            setattr(self, "dataurl",
                    self.dataurl.replace("ncep.reanalysis", "ncep.reanalysis2"))
        self.fieldlist = ["uwnd", "vwnd", "nwnd"]
        if not hasattr(self, "maxjd"):
            self.maxjd = int(pl.date2num(dtm.now())) - 3

    def setup_grid(self):
        """Setup lat-lon matrices for CCMP"""
        if not os.path.isfile(self.gridfile):
            self.retrive_file(self.dataurl, self.gridfile)
        gc = netCDF4.Dataset(self.gridfile, 'r')
        self.latvec = gc.variables['lat'][::-1].data
        self.gmt = gmtgrid.Shift(gc.variables['lon'][:].data.copy())
        self.lonvec = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    def download(self, filename):
        self.retrive_file(self.dataurl, filename)

    def load(self, field="nwnd", **kwargs):
        """Load wind field"""
        if field == "nwnd":
            uwnd = self.get_field(field="uwnd", **kwargs)
            vwnd = self.get_field(field="vwnd", **kwargs)
            setattr(self, "nwnd", np.sqrt(uwnd**2 + vwnd**2))
            return
        self._timeparams(**kwargs)

        self.filename = os.path.join(
            self.datadir, "%s.sig995.%04i.nc" % (field,self.yr))
        if not os.path.isfile(self.filename):
            self.download(self.filename)
        nc = netCDF4.Dataset(self.filename)
        jdvec = nc.variables["time"][:]/24 + self.jd0
        tpos = np.argmin(np.abs(self.jd - jdvec))
        raw = self.gmt.field(nc.variables[field][tpos,...])
        fld = raw.data
        fld[raw.mask] = np.nan
        setattr(self, field, fld[self.ijslice_flip])

    @property
    def landmask(self):
        if not hasattr(self, "_landmask"):
            gc = netCDF4.Dataset(self.gridfile)
            _landmask = gc.variables["land"][0,...].astype(bool)
            self._landmask = self.gmt.field(_landmask)[self.ijslice_flip]
        return self._landmask




























    
class __Quikscat(base.Grid):
    
    def __init__(self, **kwargs):
        super(Quikscat, self).__init__(**kwargs)
        self.fieldlist = ["uwnd", "uwnd", "nwnd"]
        
    def setup_grid(self):
        self.latvec = np.linspace(-89.875, 89.875, 720)
        self.gmt = gmtgrid.Shift(np.linspace(0, 360, 1440))
        self.lonvec = self.gmt.lonvec 
        self.llon,self.llat = np.meshgrid(self.lonvec, self.latvec)

    def download(self, filename):
        """Download missing data file from server"""
        try:
            self.retrive_file(self.dataurl, filename)
        except ftplib.error_perm:
            return False

        
    def load(self, fldname="nwnd", **kwargs):
        """Load wind field"""
        if field == "nwnd":
            uwnd = self.get_field(field="uwnd", **kwargs)
            vwnd = self.get_field(field="vwnd", **kwargs)
            setattr(self, "nwnd", np.sqrt(uwnd**2 + vwnd**2))
            return
        self._timeparams(**kwargs)


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
            
class __CORE2:

    def __init__(self,datadir = "/projData/CORE2/"):
        self.jdbase = pl.date2num(dtm(1948,1,1))+15
        self.datadir = datadir
        filename = "u_10.2005.05APR2010.nc"
        try:
            n = netCDF4.Dataset(datadir + filename)
        except:
            raise IOError("Can't open the gridfile %s" % datadir + filename)
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
            raise IOError("Can't open the gridfile %s" % datadir + filename)
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
                print(jd22-jd11,jd2-jd11,jd11,jd22)
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
