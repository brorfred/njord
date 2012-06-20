
from datetime import datetime as dtm
import math

import numpy as np
import pylab as pl
from scipy.spatial import KDTree, cKDTree
import matplotlib.cm as cm


import pycdf
from pyhdf.SD import SD,SDC

import projmaps
import gmtgrid
import figpref
from hitta import GBRY
import gmtgrid 

class NOBM:

    def __init__(self, datadir="/projData/NOBM/",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        self.i1 = 0
        self.i2 = 287
        self.j1 = 0
        self.j2 = 233
        self.y1, self.y2 = (2003,2010)
        self.m1, self.m2 = (1,12)
        self.datadir = datadir
        self.fieldpref = {'dic':'dic','ssh':'h','npp':'pp','co2':'co2'}
        lon = np.linspace(0,360,self.i2+2)[:-1]
        self.lat = np.arange(-84,72,0.667) #np.linspace(-84,72,self.j2+2)[:-1]
        self.lon,self.gr = gmtgrid.config(lon,dim=0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)

    def read_bin(self, filename):
        """ Read binary output from Poseidon """
        with open(filename) as fd:
            size = np.fromfile(fd,'<i4',count=1)[0]
            assert size == (self.i2+1) * (self.j2+1) * 4
            data = np.fromfile(fd,'<f4',count=(self.j2+1) * (self.i2+1))
            return gmtgrid.convert(
                data.reshape(self.j2+1,self.i2+1), self.gr)
            

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
        
        filename = "/mon%s%04i%02i.dat" % (self.fieldpref[fldname], yr, mn)
        self.fld = self.read_bin(self.datadir + filename)
        self.fld[self.fld>1e6] = np.nan
        self.fld[self.fld==-9999]=np.nan
        self.__dict__[fldname] = self.fld
         
    def add_ij(self):
        i1=self.i1; i2=self.i2;j1=self.j1; j2=self.j2
        self.jmat,self.imat = np.meshgrid(np.arange(self.j2-self.j1),
                                          np.arange(self.i2-self.i1))
        self.ijvec = np.vstack((np.ravel(self.imat),np.ravel(self.jmat))).T
    def add_kd(self):
        self.kd = cKDTree(list(np.vstack((np.ravel(self.llon),
                                          np.ravel(self.llat))).T))
    def ll2ij(self,lon,lat,nei=1):
        if not hasattr(self,'kd'):
            self.add_kd()
            self.add_ij()
        dist,ij = self.kd.query(list(np.vstack((lon,lat)).T),nei)
        return self.ijvec[ij][:,0],self.ijvec[ij][:,1]

    def add_landmask(self):
        g = pycdf.CDF(self.datadir + '/scb_grid.nc')
        self.landmask = g.var('mask_rho')[:]

    def add_scbij(self):
        sc = SCB()
        sci,scj = sc.ll2ij(np.ravel(self.llon),np.ravel(self.llat))
        self.sci = sci.reshape(self.i2,self.j2)
        self.scj = scj.reshape(self.i2,self.j2)

    def timeseries(self,fld,i1,j1,i2=None,j2=None):
        """Generate a time series of a field"""
        if not i2: i2 = i1+1
        if not j2: j2 = j1+1
        self.jdvec = []
        tserie = []
        for yr in np.arange(self.y1,self.y2+1):
            for mn in np.arange(self.m1,self.m2+1):
                self.load(fld,yr=yr,mn=mn)
                tserie.append(self.fld[i1:i2,j1:j2].mean())
                self.jdvec.append(pl.datestr2num('%04i%02i01' % (yr,mn)))
        return  np.array(self.jdvec), np.array(tserie)

    def movie(self,fld,k=0,jd1=733342,jd2=733342+10):
        miv = np.ma.masked_invalid

        for n,jd in enumerate(np.arange(jd1,jd2,0.25)):
            self.load(fld,jd)
            pl.clf()
            pl.pcolormesh(miv(self.__dict__[fld][0,k,:,:]))
            pl.clim(-0.1,0.1)
            pl.savefig('%s_%05i.png' % (fld,n),dpi=150)
