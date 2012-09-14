from datetime import datetime as dtm

import numpy as np
import pylab as pl
import matplotlib.cm as cm

import pycdf

import projmap
import gmtgrid
import figpref
from hitta import GBRY


class Oscar(base.Grid):
    def __init__(self, **kwargs):
        super(MODIS, self).__init__(**kwargs)

    def setup_grid(self):
        gc = pycdf.CDF(self.gridfile)
        self.lat = gc.var('latitude')[self.i1:self.i2]
        lon = gc.var('longitude')[self.j1:self.j2]
        lon[lon>360]=lon[lon>360]-360
        self.lon,self.gr = gmtgrid.config(lon, dim=0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        
    def load(self,jd=0,yr=0,mn=1,dy=1):
        """ Load Oscar fields for a given day"""
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if jd!=0:
            yr = pl.num2date(jd).year
            mn = pl.num2date(jd).month
            dy = pl.num2date(jd).day
            md = jd - pl.date2num(dtm(1992,10,05))
        elif yr!=0:
            jd = pl.date2num(yr,mn,dy)
            md  = jd - pl.date2num(dtm(1992,10,05))

        filename ="/oscar_vel%i.nc" % yr
        filenam2 ="/oscar_vel%i.nc" % (yr+1)

        nc1 = pycdf.CDF(self.datadir + filename)        
        tvec = nc1.var('time')[:]
        t1 = int(np.nonzero((tvec<=md))[0].max())
        print t1,max(tvec)
        if t1<(len(tvec)-1):
            nc2 = nc1
            t2 = t1 +1
        else:
            nc2 = pycdf.CDF(self.datadir + filenam2)                    
            t2 = 0
  
        u1 = gmtgrid.convert(nc1.var('u')[t1, 0,i1:i2,j1:j2],self.gr)
        v1 = gmtgrid.convert(nc1.var('v')[t1, 0,i1:i2,j1:j2],self.gr)
        u2 = gmtgrid.convert(nc2.var('u')[t2, 0,i1:i2,j1:j2],self.gr)
        v2 = gmtgrid.convert(nc2.var('v')[t2, 0,i1:i2,j1:j2],self.gr)
        rat = float(md-tvec[t1])/float(tvec[t2]-tvec[t1])
        self.u = u2*rat + u1*(1-rat)
        self.v = v2*rat + v1*(1-rat)
        print jd,md,t1,t2

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
        nc = pycdf.CDF(self.gridfile)
        u = nc.var('u')[:,:,:,:1081]
        msk = np.squeeze(nanmean(u,axis=0))
        msk = msk*0
        msk[np.isnan(msk)]=1
        self.landmask = gmtgrid.convert(msk,self.gr)
        
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

    def chl_mc(self):
        i1 = 0
        i2 = 480
        j1 = 0
        j2 = 1080
        mat = np.zeros( (12,i2-i1,j2-j1) )
        
        nw = pycdf.CDF('CHL_MC_AVISO.nc',pycdf.NC.WRITE|pycdf.NC.CREATE)
        nw.title = 'Monthly Climatology of MODIS Chl'
        nw.automode()
        time = nw.def_dim('Month',12)
        lat =  nw.def_dim('Latitude',i2-i1)
        lon =  nw.def_dim('Longitude',j2-j1)
        chl  = nw.def_var('Chl', pycdf.NC.FLOAT, (time, lat,lon))
        latD = nw.def_var('latitude', pycdf.NC.FLOAT, (lat,))
        lonD = nw.def_var('longitude', pycdf.NC.FLOAT, (lon,))
        timD = nw.def_var('time', pycdf.NC.FLOAT, (time,))
        
        latD[:] = self.lat[i1:i2].astype(np.float32)
        lonD[:] = self.lon[j1:j2].astype(np.float32)
        timD[:] = np.arange(1,13).astype(np.float32)


        for n,t in enumerate(np.arange(1,13)):
            fld = self.loadL3(mc=t)[i1:i2,j1:j2]
            chl[n,:,:]  = (fld).astype(np.float32)
            print t,n
        nw.close()


    def uvmat(self):
        hsmat = np.zeros ([20]+list(self.llat.shape)).astype(np.int16)
        jd1 = pl.date2num(dtm(2003,1,1))
        jd2 = pl.date2num(dtm(2009,12,31))

        vlist = np.linspace(0,1.5,21)
        for jd in np.arange(jd1,jd2+1):
            print pl.num2date(jd)
            self.load(jd=jd)
            uv = np.sqrt(self.u**2 + self.v**2)
            for n,(v1,v2) in enumerate(zip(vlist[:-1],vlist[1:])):
                msk = (uv>=v1) & (uv<v2)
                hsmat[n,msk] += 1
        return hsmat



    def movie(self):
        import matplotlib as mpl
        mpl.rcParams['axes.labelcolor'] = 'white'
        pl.close(1)
        pl.figure(1,(8,4.5),facecolor='k')
        miv = np.ma.masked_invalid
        figpref.current()
        jd0 = pl.date2num(dtm(2005,1,1))
        jd1 = pl.date2num(dtm(2005,12,31))
        mp = projmaps.Projmap('glob')
        x,y = mp(self.llon,self.llat)
        for t in np.arange(jd0,jd1):
            print pl.num2date(t)
            self.load(t)
        
            pl.clf()
            pl.subplot(111,axisbg='k')
            mp.pcolormesh(x,y,
                          miv(np.sqrt(self.u**2 +self.v**2)),
                          cmap=cm.gist_heat)
            pl.clim(0,1.5)
            mp.nice()
            pl.title('%04i-%02i-%02i' % (pl.num2date(t).year,
                                         pl.num2date(t).month,
                                         pl.num2date(t).day),
                     color='w')
            pl.savefig('/Users/bror/oscar/norm/%03i.png' % t,
                       bbox_inches='tight',facecolor='k',dpi=150)
        
