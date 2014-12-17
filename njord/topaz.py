import os
import glob
import datetime
from datetime import datetime as dtm
import re

import numpy as np
import pylab as pl
import random
import matplotlib as mpl
from matplotlib.mlab import griddata
import matplotlib.cm as cm
from scipy.stats import nanmean
from scipy.io import netcdf_file

from netCDF4 import Dataset

import base
import winds
import gmtgrid
import grid
import bln
from pysea import sw, wind2pv
from reuerflux import pv2fl, pv2wpv
lsd = 3600*24*1000

def nodimb(a,b): return (a-b) / b


class Tpz(base.Grid):
    def __init__(self, **kwargs):
        """Init the class and setup the grid. Paths are defined here."""
        super(Tpz, self).__init__(**kwargs)
        self.lat = self.lat[self.j1:self.j2].copy()
        self.lon = self.lon[self.i1:self.i2].copy()
        self.create_pardict(lnmsk='', dtmsk='')
        self.keypref = '1d_1yr_'

        self._setup_time()
        self._setup_fields()
        self.gcm = 'tpz'


    def setup_grid(self):
        """Setup necessary variables for grid """
        g = netcdf_file(self.gridfile)

        self.gmt = gmtgrid.Shift(g.variables['gridlon_t'][:].copy())
        self.lon = (self.gmt.lonvec).astype(np.float32).copy()
        self.lat   = g.variables['gridlat_t'][:].copy()
        self.depth = self.gmt.field(g.variables['depth_t'][:])
        self.dxt   = self.gmt.field(g.variables['dxt'][:])
        self.dyt   = self.gmt.field(g.variables['dyt'][:])
        self.dxu   = self.gmt.field(g.variables['dxu'][:])
        self.dyu   = self.gmt.field(g.variables['dyu'][:])
        #self.dzt   = self.gmt.field(g.variables['dz_t'][:])
        #self.zlev  = self.gmt.field(g.variables['dz_t'][:])
        #if not hasattr(self, 'k1'): self.k1 = 0
        #if not hasattr(self, 'k2'): self.k2 = len(self.zlev)
        #self.dz    =  g.variables('dz_t')[self.k1:self.k2,25,30]
        #self.vol = ( self.dxt[np.newaxis,...] * 
        #             self.dyt[np.newaxis,...]*self.dzt)
        #self.vol[self.vol<0] = np.nan
        #self.zlev = [0,] + np.cumsum(self.dz)-10
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.area = self.dxt * self.dyt
        #for v in ['dzt','llon','llat']:
        #    self.__dict__[v] = self.__dict__[v].astype(np.float32)


    def _setup_fields(self):
        self.pa = {'nwnd':'wind','mldp':'mld','o2st':'o2_saturation',
                   'o2ct':'o2', 'arst':'ar_sat_o2satar', 'ncpo':'jo2',
                   'arct':'ar', 'temp':'temp', 'salt':'salt',
                   'hblt':'hblt','no3c':'no3','fedc':'fed',
                   'poco':'nphyto_tot',
                   'o2fl':'sfc_flux_o2','arfl':'sfc_flux_ar',
                   'cofl':'sfc_flux_co2',
                   'neo2':'neutral_o2','dfo2':'o2_zflux_diff',
                   'ssht':'eta_t','chlo':'chl',
                   'uvel':'u','vvel':'v','sshu':'dhu'}

        self.slope={}
        for a in self.pa:
            self.slope[a]=1
        self.slope['ncpo'] = lsd
        self.slope['o2fl'] = -lsd
        self.slope['cofl'] = -lsd
        self.slope['neo2'] = lsd
        self.slope['dfo2'] = lsd
        self.slope['o2ct'] = 1000 
        self.slope['o2st'] = 1000
        self.slope['arct'] = 1000
        self.slope['arst'] = 1000
        self.slope['poco'] = 6.625*1000
        self.slope['no3c'] = 1000
        self.slope['fedc'] = 1e6
        self.slope['uvel'] = 3600*24
        self.slope['vvel'] = 3600*24
        self.fields = []


    def _setup_time(self):

        if not hasattr(self, 't1'): self.t1 = 0
        if not hasattr(self, 't2'): self.t2 = None
        def load(keysuff):
            print self.vd[self.keypref + keysuff][2]
            n = Dataset(self.vd[self.keypref + keysuff][2])
            return n.variables[self.vd[self.keypref + keysuff][0]][:]
        self.tvec = load('time')[self.t1:self.t2] + 623542
        self.tops = self.tvec

    def condload(self, parlist):
        for par in (parlist):
            if not hasattr(self,par):
                self.load(par)

    def load(self, par, surf=False,mldint=False,fullint=False,noneg=False):
        self.fields.append(par)

        def scale(fld):
            return (fld*self.slope[par]).astype(np.float32)
        if par == 'dens':
            self.load('temp',surf,mldint,fullint)
            self.load('salt',surf,mldint,fullint)
            self.dens = sw.dens0(self.salt, self.temp)
            return
        elif par == 'o2ar':
            self.condload( ('o2ct','o2st','arct','arst') )
            self.o2ar =  ( (self.o2ct/self.o2st) /
                           (self.arct/self.arst) - 1 ) * 100
            return
        elif par == 'd003':
            self.condload( ('mldp','dens') )
            mld = self.mldp.copy() * 0
            for k in np.arange(len(self.zlev)):
                mld[(self.dens[:,k,...]-self.dens[:,2,...])<0.03] = k
            for k in np.arange(len(self.zlev),0,-1)-1:
                mld[mld==k] = self.zlev[k]-self.dz[k]/2.0
            self.d003 = mld
            return
        elif par == 'depth':
            return
        elif par == 'mldk':
            self.condload( ('mldp',) )
            self.mldk = np.ceil(np.interp(self.mldp, 
                                          self.zlev, np.arange(50)) )
            return
  
        
        elif par == 'wwfl':
            self.condload( ('o2ar','dens','nwnd','mldp') )
            del self.o2ct, self.arct, self.arst
            self.wwfl = pv2fl(self.o2ar, self.temp, self.mldp, self.nwnd,
                              wtlen=60,o2st=self.o2st, dens=self.dens)
            return
        elif par == 'crwwfl':
            self.fields[-1] = 'wwfl'
            tm = len(self.tvec)
            tvec = self.tvec
            self.tvec = np.arange(self.tvec.min()-61,self.tvec.min())
            self.load('crwnd')
            self.tvec = tvec
            self.condload( ('nwnd','o2ar','dens','mldp') )
            self.nwnd = np.concatenate((self.crwnd,self.nwnd),axis=0)
            fillmat = np.zeros( (self.crwnd.shape[0],) + self.o2ar.shape[1:] )
            for v in ['o2ar','temp','o2st','dens','mldp']:
                self.__dict__[v] = np.concatenate(
                    (self.__dict__[v][tm-60:,...],self.__dict__[v]), axis=0)
            try:
                del self.o2ct, self.arct, self.arst,self.crwnd
            except AttributeError:
                pass
            print self.nwnd.shape,self.o2ar.shape
            self.wwfl = pv2fl(self.o2ar, self.temp, self.mldp, self.nwnd,
                              wtlen=60,o2st=self.o2st, dens=self.dens)
            for v in  ['o2ar','temp','mldp','o2st','dens','wwfl','nwnd']:
                self.__dict__[v] = self.__dict__[v][60:,...] 
            return
        elif par == 'nn10':
            if not hasattr(self,'ncpm'):
                self.load('ncpo',mldint=True)
            self.nn10 =self.ncpm.copy() * np.nan
            for t in np.arange(10,self.nn10.shape[0]):
                self.nn10[t,...]=np.mean(self.ncpm[t-10:t,...],axis=0)
            return
        elif par == 'nnlg':
            if not hasattr(self,'ncpm'):
                self.load('ncpo',mldint=True)
            self.condload( ('mldp','wpv','temp') )
            self.nnlg = self.ncpm.copy()
            self.nnlg[self.nnlg>1e6] = 0
            self.nnlg[self.nnlg<-1e6] = 0
            self.cnt  = self.ncpm[:] * 0
            wlg = np.floor((self.mldp/self.wpv)/np.log(2))
            wlg[wlg>60] = 60
            for t in np.arange(60,self.ncpm.shape[0]):
                for b in np.arange(60):
                    msk = wlg[t,:,:]>b
                    self.nnlg[t,msk] = self.nnlg[t,msk] + self.ncpm[t-b,msk]
                    self.cnt[t,msk]  = self.cnt[t,msk] + 1
            self.nnlg = self.nnlg/self.cnt
            return
        
        elif par == 'nrm10':
            self.condload( ('crwwfl','nn10') )
            self.nrm10 = nodimb(self.wwfl,self.nn10)
            return
        elif par == 'wpv':
            self.condload( ('wwfl',) )
            self.wpv = pv2wpv(self.temp,self.mldp,self.nwnd,dens=self.dens)
            return
        elif par == 'crwpv':
            tm = len(self.tvec)
            tvec = self.tvec
            self.tvec = np.arange(self.tvec.min()-61,self.tvec.min())
            self.load('crwnd')
            self.tvec = tvec
            self.condload( ('nwnd','o2ar','dens','mldp') )
            self.nwnd = np.concatenate((self.crwnd,self.nwnd),axis=0)
            fillmat = np.zeros( (self.crwnd.shape[0],) + self.o2ar.shape[1:] )
            for v in ['temp','dens','mldp']:
                print v, " before: ",self.__dict__[v].shape
                self.__dict__[v] = np.concatenate(
                    (self.__dict__[v][tm-60:,...],self.__dict__[v]), axis=0)
                print v, " after: ",self.__dict__[v].shape
            self.wpv = pv2wpv(self.temp,self.mldp,self.nwnd,dens=self.dens)
            for v in  ['temp','mldp','dens','wpv','nwnd']:
                self.__dict__[v] = self.__dict__[v][60:,...] 
            return
        
        elif par == 'crwnd':
            self.fields[-1] = 'nwnd'
            cr = winds.CORE2()
            cr.load(self.tvec.min(), self.tvec.max())
            self.crwnd = np.zeros(self.tvec.shape + self.llon.shape)
            for tpos in np.arange(len(self.tvec)):
                try:
                    self.crwnd[tpos,:,:] = self.reproject(cr, cr.nwnd[tpos,:,:])
                except IndexError:
                    self.crwnd[tpos,:,:] = self.reproject(cr, cr.nwnd[-1,:,:])
            return
        
        key = self.keypref + self.pa[par]
        self.par=par; vd=self.vd;
        t1 = self.t1; t2 = self.t2
        i1 = self.i1; i2 = self.i2
        j1 = self.j1; j2 = self.j2
        k1 = self.k1; k2 = self.k2        
        if not t1: t1 = 0
        if not t2: t2 = self.tvec.shape[0]
        n = netcdf_file(vd[key][2])

        if par == 'nwnd':
            fld = self.gmt.field(scale(n.variables[vd[key][0]][:])[t1:t2,...])
            fld = fld#[:,i1:i2,j1:j2]
            fld[fld <-1e9] = np.nan
            fld[fld > 1e9] = np.nan
            self.__dict__[par] = fld.astype(np.float32)
            return  
        if  'zt_ocean' in vd[key][1].dimensions and mldint:
            """ Integrate down to MLD """
            na = np.newaxis
            exists = False
            if not hasattr(self,par): exists = True
            self.condload( (par,'mldp',) )
            fld = grid.mldint(self.__dict__[par],self.mldp,self.dz)
            if exists: del self.__dict__[par]
            par = par[:-1] + 'm'
            self.fields[-1] = par

        elif  'zt_ocean' in vd[key][1].dimensions and fullint:
            """ Integrate over each column """
            na = np.newaxis
            exists = False
            if not hasattr(self,par): exists = True
            self.condload( (par,'mldp','ssht') )
            fld = self.__dict__[par] * self.dz[na,:,na,na]
            if noneg: fld[fld<0] = 0
            fld[:,0,:,:] = fld[:,0,:,:] * (self.ssht+self.dz[0]) / self.dz[0]
            fld = np.nansum(fld,axis=1)
            if exists: del self.__dict__[par]
            par = par[:-1] + 'f'
            self.fields[-1] = par

        elif 'zt_ocean' in vd[key][1].dimensions and surf:
            fld = self.gmt.field(scale(n.variables[vd[key][0]][t1:t2,0,    j1:j2,:]))
        elif 'zt_ocean' in vd[key][1].dimensions:
            fld = self.gmt.field(scale(n.variables[vd[key][0]][t1:t2,k1:k2,j1:j2,:]))
        else:
            fld = self.gmt.field(scale(n.variables[vd[key][0]][t1:t2,      j1:j2,:]))
        print par + ": ", fld.shape
        fld[fld <-1e9] = np.nan
        fld[fld > 1e9] = np.nan
        self.__dict__[par] = fld[...,self.i1:self.i2]

    def t2iso(t):
        baseiso = mpl.dates.date2num(datetime.datetime(1948,1,1))
        return baseiso + t - 1

    def create_pardict(self,lnmsk='',dtmsk=''):
        vardict = {}
        for f in glob.glob(self.datadir + "/*/*.nc"):
            if not "month" in f:
                n = Dataset(f)
                for k in n.variables.keys():
                    ln,dt = f.split('/')[-2].split('_')
                    if  ( (ln==lnmsk or lnmsk=='') and
                          (dt==dtmsk or dtmsk=='') ): 
                        fld_index = dt + "_" + ln + "_" + k
                        vardict[fld_index] = ( (k,) + (n.variables[k],)
                                               + (f,) + (ln,)+(dt,) )
        self.vd = vardict
    
    def attr(self,key):
        n = netcdf_file(self.vd[key][2])
        at = n.variables[self.vd[key][0]]
        return at.attributes()

    def scl(self,a):
        f = np.load('result/windlatlon.npz')
        llat = np.ravel(f['llat'])
        llon = np.ravel(f['llon'])
        def imr(a):
            return griddata(llat, llon, np.ravel(a),
                            self.llat, self.llon, interp='linear' )
        if a.ndim > 2 :
            b = np.zeros((a.shape[-3],len(self.lat),len(self.lon)))
            for t in np.arange(a.shape[-3]):
                b[t,:,:]=imr(a[t,:,:])
        else:
            b = imr(a)
        return b

    def add_year(self):
        if not hasattr(self,'year_added'): 
            self.year_added=[]
            self.tvec = np.concatenate((self.tvec,self.tvec+365),axis=0)
        for fld in self.fields:
            if not fld in self.year_added and fld in self.__dict__:
                self.__dict__[fld] = np.concatenate( 
                    (self.__dict__[fld],self.__dict__[fld]) ,axis=0)
                self.year_added.append(fld)

    def create_months(self):
        self.months = np.array([ dt.month for dt in pl.num2date(self.tvec)])
        
    def get_corewind(self):
        xi,yi = np.meshgrid(np.arange(360),np.arange(200))
        x,y = np.meshgrid(np.arange(94),np.arange(192))

    def add_landmask(self):
        self.landmask = self.depth == 0

    def add_seaslen(self,dfn='ncpm'):
        if dfn == 'o2ar':
            self.condload( ('o2ar',) )
            if self.o2ar.ndim == 3:
                msk = self.o2ar[:,:,:].copy()
            else:
                msk = self.o2ar[:,0,:,:].copy()
            msk[msk>0] = 1
            msk[msk<0] = 0
            msk[np.isnan(msk)] = 0
            seaslen = np.sum(msk,axis=0) * 365 / len(self.tvec)
            self.seaslen = msk
            self.seaslen[:] = seaslen
        else:
            if not hasattr(self,'ncpm'):
                self.load('ncpo',mldint=True)
            msk = self.ncpm[:,:,:].copy()
            msk[msk<=0] = 0
            msk[msk>0] = 1
            msk[np.isnan(msk)] = 0
            seaslen = np.sum(msk,axis=0) * 365 / len(self.tvec)
            self.seaslen = msk
            self.seaslen[:] = seaslen

