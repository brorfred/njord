
import pylab as pl
import numpy as np
from scipy.io import netcdf_file

import base



class SCB(base.Grid):
    """ Manipulate data from the official jpl SCB runs """
    def __init__(self, datadir="/projData/jplSCB/ROMS/",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        self.i1 = 0
        self.i2 = 111
        self.j1 = 0
        self.j2 = 211
        self.datadir = datadir
        g = netcdf_file(datadir + '/scb_das_grid.nc', SDC.READ)
        self.lat = q.variables['lat'][:]
        self.lon = g.variables['lon'][:]-360
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)

    def load(self,fldname,jd=0,yr=0,mn=1,dy=1,hr=3):
        """ Load Cali Current fields for a given day"""
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if jd!=0:
            yr = pl.num2date(jd).year
            mn = pl.num2date(jd).month
            dy = pl.num2date(jd).day
            hr = pl.num2date(jd+0.125).hour
            mi = pl.num2date(jd).minute
            se = pl.num2date(jd).second      
        elif yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy,hr))
        yd = jd - pl.date2num(dtm(yr,1,1)) + 1
        filename = "/scb_das_%04i%02i%02i%02i.nc" % (yr,mn,dy,hr)
        print filename
        nc = pycdf.CDF(self.datadir + filename)        
        fld =  nc.var(fldname)[:]
        fld[fld==-9999]=np.nan
        self.__dict__[fldname] = fld
         
    def add_landmask(self):
        nc = pycdf.CDF(self.datadir + '/scb_das_grid.nc')
        msk = nc.var('temp')[0,0,...]
        msk[msk!=-9999] = 1
        msk[msk==-9999] = 0
        self.landmask = msk

 
class NOW(base.Grid):
    """ Manipulate data from the NOWCAST jpl SCB runs """
    def __init__(self, **kwargs):
        super(NOW, self).__init__()
        self.add_mp()

    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile)
        self.llat = g.variables['lat_rho'][:].copy()
        self.llon = g.variables['lon_rho'][:].copy()
        self.depth = g.variables['h'][:].copy()
        self.Cs_r = np.array([-0.882485522505154, -0.778777844867132,
                              -0.687254423585503, -0.606483342183883,
                              -0.535200908367393, -0.472291883107274,
                              -0.416772032329648, -0.367772728223386,
                              -0.324527359249072, -0.286359336228826,
                              -0.252671506867986, -0.222936813095075,
                              -0.196690045050832, -0.173520562714503,
                              -0.153065871294677, -0.135005949869352,
                              -0.11905824454479,  -0.104973247799366,
                              -0.0925305948496471, -0.0815356159649889,
                              -0.0718162907903607, -0.0632205570267136,
                              -0.0556139313622304, -0.0488774054330831,
                              -0.042905583895255, -0.0376050354769002,
                              -0.0328928312128658, -0.0286952469915389,
                              -0.0249466101148999, -0.021588271825806,
                              -0.0185676897273263, -0.0158376057382559,
                              -0.0133553067236276, -0.0110819562325242,
                              -0.00898198688799211, -0.00702254392277909,
                              -0.00517297115481568, -0.0034043313603439,
                              -0.00168895354075999, 0.])

    def load(self,fldname, **kwargs):
        """ Load Cali Current fields for a given day"""
        if self._calcload(fldname, **kwargs): return
        self._timeparams(**kwargs)
        
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        filename = ("/%04i%02i%02i%02i_da.nc" % 
                    (self.yr,self.mn,self.dy,self.hr) )
        print self.datadir + filename
        nc = netcdf_file(self.datadir + filename)        
        fld =  np.squeeze(nc.variables[fldname][:]).copy()
        fld[fld==-9999]=np.nan
        self.__dict__[fldname] = fld
        self.ssh =  np.squeeze(nc.variables['zeta'][:])
        self.zlev = ((self.depth + self.ssh)[np.newaxis,:,:] *
                     self.Cs_r[:,np.newaxis,np.newaxis])

    @property
    def landmask(self):
        """ Generate a landmask for t-positions"""
        return self.depth <= 10
