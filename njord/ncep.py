from datetime import datetime as dtm

import numpy as np
import pylab as pl
from scipy.spatial import KDTree, cKDTree
import matplotlib.cm as cm


import pycdf

import projmaps
import gmtgrid
import figpref
from hitta import GBRY

class Reanalysis:

    def __init__(self, datadir="/projData/ncep",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        self.i1 = 0
        self.i2 = 73
        self.j1 = 0
        self.j2 = 144
        self.datadir = datadir
        gc = pycdf.CDF(self.datadir + '/land.nc')
        self.lat = gc.var('lat')[self.i1:self.i2]
        lon = gc.var('lon')[self.j1:self.j2]
        lon[lon>360]=lon[lon>360]-360
        self.lon,self.gr = gmtgrid.config(lon, dim=0)
        self.llon,self.llat = np.meshgrid(self.lon,self.lat)
        
    def load(self,field, jd=0,yr=0,mn=1,dy=1):
        """ Load NCEP reanalysis fields for a given day"""
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
        return self.ijvec[ij-1][:,0],self.ijvec[ij-1][:,1]

    def add_landmask(self):
        """Load and add landmask"""
        nc = pycdf.CDF(self.datadir + '/land.nc')
        self.landmask = gmtgrid.convert( nc.var('land')[:],self.gr)
