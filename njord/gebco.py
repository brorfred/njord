from datetime import datetime as dtm

import numpy as np
import pylab as pl
import matplotlib.cm as cm
from scipy.io import netcdf_file

import gmtgrid
#import figpref
#from hitta import GBRY

import base

class Thirty(base.Grid):
    def __init__(self, **kwargs):
        super(Thirty, self).__init__(**kwargs)

    def setup_grid(self):
        gc = netcdf_file(self.gridfile)
        dlon,dlat = gc.variables['spacing'][:]
        lon1,lon2 = gc.variables['x_range'][:]
        lat1,lat2 = gc.variables['y_range'][:]
        self.lon = np.arange(lon1, lon2, dlon)
        self.lat = np.arange(lat1, lat2, dlat)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        
    def load(self,**kwargs):
        """ Load GEBCO bathymetry """
        z = nc.variables['z']
        self._timeparams(**kwargs)
        md  = self.jd - pl.datestr2num('1992-10-05')
        
    def init_nasa(self):
        import nasa
        self.ns = nasa.MODIS(res='9km')

    def loadL3(self,jd=732281.0, mc=0, fld='chl'):
        if not hasattr(self,'gcmi'):
            self.create_satijvecs()
        if mc != 0:
            self.ns.load(fld,fldtype="MC",mc=mc)
        else:
            yr = pl.num2date(jd).year
            yd = jd - pl.date2num(dtm(pl.num2date(jd).year,1,1)) + 1
            self.ns.load(fld,yr=yr,yd=yd)
        cnt = self.llat[:] * 0
        fld = self.llat[:] * 0
        msk = ~np.isnan(np.ravel(self.ns.chl))
        for gj,gi,sj,si in np.vstack((self.satj[msk],self.sati[msk],
                                      self.gcmj[msk],self.gcmi[msk])).T: 
            fld[gj,gi] += self.ns.chl[sj,si]
            cnt[gj,gi] += 1
        fld = fld/cnt
        return fld

    def add_landmask(self):
        from scipy.stats import nanmean
        self.load()
        self.landmask = np.isnan(self.u)
        
    def dchl_dt_time(self):
        i1 = 0
        i2 = 230
        j1 = 275
        j2 = 600
        t1 = pl.date2num(dtm(2009,1,1))
        t2 = pl.date2num(dtm(2009,12,31))
        mat = np.zeros( (t2-t1,i2-i1,j2-j1) )

        nw = pycdf.CDF('dchl2009.nc',pycdf.NC.WRITE|pycdf.NC.CREATE)
        nw.title = 'DChl/Dt'
        nw.automode()
        time = nw.def_dim('time',t2-t1)
        lat =  nw.def_dim('Latitude',i2-i1)
        lon =  nw.def_dim('Longitude',j2-j1)
        dchl = nw.def_var('DChlDt', pycdf.NC.FLOAT, (time, lat,lon))
        chl  = nw.def_var('Chl', pycdf.NC.FLOAT, (time, lat,lon))
        u    = nw.def_var('u', pycdf.NC.FLOAT, (time, lat,lon))
        v    = nw.def_var('v', pycdf.NC.FLOAT, (time, lat,lon))
        latD = nw.def_var('latitude', pycdf.NC.FLOAT, (lat,))
        lonD = nw.def_var('longitude', pycdf.NC.FLOAT, (lon,))
        timD = nw.def_var('time', pycdf.NC.FLOAT, (time,))
        
        latD[:] = self.lat[i1:i2].astype(np.float32)
        lonD[:] = self.lon[j1:j2].astype(np.float32)
        timD[:] = np.arange(t1,t2).astype(np.float32)

        fld1 = self.loadL3(t1)[i1:i2,j1:j2]
        for n,t in enumerate(np.arange(t1+1,t2)):
            fld2 = self.loadL3(t)[i1:i2,j1:j2]
            self.load(t-1)            
            dchl[n,:,:] = (fld2-fld1).astype(np.float32)
            chl[n,:,:]  = (fld1).astype(np.float32)
            u[n,:,:]    = (self.u).astype(np.float32)[i1:i2,j1:j2]
            v[n,:,:]    = (self.v).astype(np.float32)[i1:i2,j1:j2]
            fld1 = fld2
            print t,n
        nw.close()

