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

import pycdf

import gmtgrid
import grid
import bln
from pysea import sw, wind2pv
from reuerflux import pv2fl, pv2wpv
lsd = 3600*24*1000

def nodimb(a,b): return (a-b) / b


class Tpz:
    def __init__(self, datadir='/projData/TOPAZ/',
                 t1=None, t2=None, k1=None, k2=None,
                 i1=None, i2=None, j1=None, j2=None,area=None,
                 ln='', dt=''):
        """Init the class and setup the grid. Paths are defined here."""
        self.datadir = datadir
        self.t1 = t1; self.t2 = t2
        self.i1 = i1; self.i2 = i2
        self.j1 = j1; self.j2 = j2
        self.k1 = k1; self.k2 = k2
        if area:
            self.i1 = area(0); self.i2 = area(1)
            self.j1 = area(2); self.j2 = area(3)
            self.k1 = area(4); self.k2 = area(5)        
        self.create_pardict(ln, dt)
        self.keypref = '1d_1yr_'

        vd = self.vd
        def load(keysuff):
            n = pycdf.CDF(vd[self.keypref + keysuff][5])
            return n.var(vd[self.keypref + keysuff][0])[:]

        self.tvec = load('time')[t1:t2] + 623542
        self.tops = self.tvec

        #Grid definitions
        g = pycdf.CDF(datadir + 'grid_spec_cm2_ocean_tripolar.nc')

        lon = g.var('gridlon_t')[:]
        self.lon,self.gr = gmtgrid.config(lon,dim=0)
        self.lon = self.lon[j1:j2]
        self.lat = g.var('gridlat_t')[i1:i2]
        self.dxt = gmtgrid.convert(g.var('dxt')[:], self.gr)[i1:i2,j1:j2]
        self.dyt = gmtgrid.convert(g.var('dyt')[:], self.gr)[i1:i2,j1:j2]
        self.dxu = gmtgrid.convert(g.var('dxu')[:], self.gr)[i1:i2,j1:j2]
        self.dyu = gmtgrid.convert(g.var('dyu')[:], self.gr)[i1:i2,j1:j2]
        self.dzt = gmtgrid.convert(g.var('dz_t')[:], self.gr)[k1:k2,i1:i2,j1:j2]
        self.zlev = gmtgrid.convert(g.var('dz_t')[:], self.gr)[k1:k2]
        self.dz =  g.var('dz_t')[k1:k2,25,30]
        self.vol = ( self.dxt[np.newaxis,...] * 
                     self.dyt[np.newaxis,...]*self.dzt)
        self.vol[self.vol<0] = np.nan
        self.zlev = [0,] + np.cumsum(self.dz)-10
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        self.area = self.dxt * self.dyt

        for v in ['dzt','llon','llat']:
            self.__dict__[v] = self.__dict__[v].astype(np.float32)

        #Loadable fields
        self.pa = {'nwnd':'wind','mldp':'mld','o2st':'o2_saturation',
                   'o2ct':'o2', 'arst':'ar_sat_o2satar', 'ncpo':'jo2',
                   'arct':'ar', 'temp':'temp', 'salt':'salt',
                   'hblt':'hblt','no3c':'no3','fedc':'fed',
                   'poco':'nphyto_tot',
                   'o2fl':'sfc_flux_o2','arfl':'sfc_flux_ar',
                   'neo2':'neutral_o2','dfo2':'o2_zflux_diff',
                   'ssht':'eta_t','chlo':'chl',
                   'uvel':'u','vvel':'v','sshu':'dhu'}

        self.slope={}
        for a in self.pa:
            self.slope[a]=1
        self.slope['ncpo'] = lsd
        self.slope['o2fl'] = -lsd
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

        self.gcm = 'tpz'

        def lim(val1,val2):
            if not val1: return val2
            return val1
        self.t1 = lim(t1,0); self.t2 = lim(t2,len(self.tvec))
        self.i1 = lim(i1,0); self.i2 = lim(i2,self.llat.shape[0])
        self.j1 = lim(j1,0); self.j2 = lim(j2,self.llat.shape[1])
        self.k1 = lim(k1,0); self.k2 = lim(k2,len(self.zlev))

        self.fields = []

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
        elif par == 'mldk':
            self.condload( ('mldp',) )
            self.mldk = np.ceil(np.interp(self.mldp, 
                                          self.zlev, np.arange(50)) )
            return
            """
            if d003:
                if not hasattr(self,'d003'):
                    self.load('d003')
                mld = self.d003
            else:
                self.condload( ('o2ct','o2st','arct','arst') )

            mld = self.mldp
            zl = np.cumsum(self.dz)
            self.mldk = np.interp(mld, zl, np.arange(50))
            """
        
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
            
            #wlg = np.floor((self.mldp/
            #                wind2pv(self.temp[:,0,:,:],self.nwnd))/np.log(2))
            wlg = np.floor((self.mldp/self.wpv)/np.log(2))
            wlg[wlg>60] = 60
            for t in np.arange(60,self.ncpm.shape[0]):
                for b in np.arange(60):
                    msk = wlg[t,:,:]>b
                    self.nnlg[t,msk] = self.nnlg[t,msk] + self.ncpm[t-b,msk]
                    self.cnt[t,msk]  = self.cnt[t,msk] + 1
            self.nnlg = self.nnlg/self.cnt
            #del self.cnt
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
            import winds
            cr = winds.core2()
            a = cr.load(self.tvec.min(),self.tvec.max())
            llat = np.ravel(cr.llat)
            llon = np.ravel(cr.llon)
            def imr(a):
                return griddata(llat, llon, np.ravel(a),
                                self.llat, self.llon)
            if a.ndim > 2 :
                b = np.zeros((a.shape[-3],len(self.lat),len(self.lon)))
                for t in np.arange(a.shape[-3]):
                    b[t,:,:]=imr(a[t,:,:])
            else:
                b = imr(a)
            self.crwnd = b
            return

        key = self.keypref + self.pa[par]
        self.par=par; vd=self.vd;
        t1 = self.t1; t2 = self.t2
        i1 = self.i1; i2 = self.i2
        j1 = self.j1; j2 = self.j2
        k1 = self.k1; k2 = self.k2        
        if not t1: t1 = 0
        if not t2: t2 = self.tvec.shape[0]
        n = pycdf.CDF(vd[key][5])

        if par == 'nwnd':
            fld = self.utcgrid(scale(n.var(vd[key][0])[t1:t2,...]))
            fld = fld#[:,i1:i2,j1:j2]
            fld[fld <-1e9] = np.nan
            fld[fld > 1e9] = np.nan
            self.__dict__[par] = fld.astype(np.float32)
            return  
        if  'zt_ocean' in vd[key][1] and mldint:
            """ Integrate down to MLD """
            na = np.newaxis
            exists = False
            if not hasattr(self,par): exists = True
            self.condload( (par,'mldp',) )
            fld = grid.mldint(self.__dict__[par],self.mldp,self.dz)
            if exists: del self.__dict__[par]
            par = par[:-1] + 'm'
            self.fields[-1] = par

        elif  'zt_ocean' in vd[key][1] and fullint:
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

        elif 'zt_ocean' in vd[key][1] and surf:
            fld = self.utcgrid(
                scale(n.var(vd[key][0])[t1:t2,0,i1:i2,j1:j2]))
        elif 'zt_ocean' in vd[key][1]:
            fld = self.utcgrid(
                scale(n.var(vd[key][0])[t1:t2,k1:k2,i1:i2,j1:j2]))
        else:
            fld = self.utcgrid(
                scale(n.var(vd[key][0])[t1:t2,i1:i2,j1:j2]))
        print par + ": ", fld.shape
        fld[fld <-1e9] = np.nan
        fld[fld > 1e9] = np.nan
        self.__dict__[par] = fld

    def t2iso(t):
        baseiso = mpl.dates.date2num(datetime.datetime(1948,1,1))
        return baseiso + t - 1

    def utcgrid(self, fld):
        if hasattr(self,'par') and self.par == "nwnd":
            fld = self.scl(fld)
        elif fld.shape[-1]!=360:
            pass
        else:
            wstfld=fld[...,:100].copy()
            fld[...,:260] = fld[...,100:]
            fld[...,260:] = wstfld
            del wstfld
        return fld
    
    def create_pardict(self,lnmsk='',dtmsk=''):
        vardict = {}
        for f in glob.glob(self.datadir + "/*/*.nc"):
            if not "month" in f:
                n = pycdf.CDF(f)
                for k in n.variables().keys():
                    ln,dt = f.split('/')[-2].split('_')
                    if  ( (ln==lnmsk or lnmsk=='') and
                          (dt==dtmsk or dtmsk=='') ): 
                        fld_index = dt + "_" + ln + "_" + k
                        vardict[fld_index] = ( (k,) + n.variables()[k]
                                               + (f,) + (ln,)+(dt,) )
        self.vd = vardict
    
    def attr(self,key):
        n = pycdf.CDF(self.vd[key][5])
        at = n.var(self.vd[key][0])
        return at.attributes()

    def scl(self,a):
        f = np.load('result/windlatlon.npz')
        llat = np.ravel(f['llat'])
        llon = np.ravel(f['llon'])
        def imr(a):
            return griddata(llat, llon, np.ravel(a),
                            self.llat, self.llon)
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


    def add_ij(self):
        self.jmat,self.imat = np.meshgrid(np.arange(self.j2-self.j1),
                                          np.arange(self.i2-self.i1))
        self.ijvec = np.vstack((np.ravel(self.imat),np.ravel(self.jmat))).T
    def add_kd(self,mask=None):
        from scipy.spatial import KDTree, cKDTree

        latvec = np.ravel(self.llat)
        lonvec = np.ravel(self.llon)
        if not mask is None: 
            latvec = latvec[~np.isnan(np.ravel(mask))]
            lonvec = lonvec[~np.isnan(np.ravel(mask))]
            self.add_ij()
            self.ijvec = self.ijvec[~np.isnan(np.ravel(mask))]
        self.kd = cKDTree(list(np.vstack((lonvec,latvec)).T))
    def ll2ij(self,lon,lat,nei=1):
        if not hasattr(self,'imat'):
            self.add_kd()
            self.add_ij()
        dist,ij = self.kd.query(list(np.vstack((lon,lat)).T),nei)
        return self.ijvec[ij-1][:,0],self.ijvec[ij-1][:,1]
