from datetime import datetime as dtm

import numpy as np
import pylab as pl
from scipy.spatial import KDTree, cKDTree
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize

import pycdf
from pyhdf.SD import SD,SDC

import projmaps
import gmtgrid
import figpref
from hitta import GBRY

class Cal:

    def __init__(self, datadir="/projData/sat/mati/cal/",ijarea=[],
                 lat1=None,lat2=None,lon1=None,lon2=None):
        self.i1 = 0
        self.i2 = 3405
        self.j1 = 0
        self.j2 = 3840
        self.datadir = datadir
        g = SD(datadir + 'cal_aco_3840_Latitude_Longitude.hdf', SDC.READ)
        self.llat = g.select('Latitude')[:]
        self.llon = g.select('Longitude')[:]
        
    def load(self,fldname,jd=0,yr=0,mn=1,dy=1):
        """ Load Cali Current fields for a given day"""
        i1=self.i1; i2=self.i2; j1=self.j1; j2=self.j2
        if jd!=0:
            yr = pl.num2date(jd).year
            mn = pl.num2date(jd).month
            dy = pl.num2date(jd).day
        elif yr!=0:
            jd = pl.date2num(dtm(yr,mn,dy))
        md = jd - pl.date2num(dtm(1992,10,05))
        yd = jd - pl.date2num(dtm(yr,1,1)) + 1
        
        if fldname == 'chl':
            filename = "/C%04i%03i_chl_mapped.hdf" % (yr,yd)
            #ncfieldname = 'chl_%04i_%03i' % (yr,yd)
            def scale(PV): return 10**(PV*0.015-2)
        elif fldname == 'sst':
            filename = "/M%04i%03i_sst_mapped.hdf" % (yr,yd)
            #ncfieldname = 'sst_%04i_%03i' % (yr,yd)            
            def scale(PV): return PV*0.15000001-3
        h = SD(self.datadir + filename,SDC.READ)        
        ncfieldname = h.datasets().keys()[0]
        fld =  h.select(ncfieldname)
        attr = fld.attributes()
        PV = fld[:].astype(np.float)
        PV[PV<0] = PV[PV<0]+256
        PV[PV==0]   = np.nan
        PV[PV==255] = np.nan
        self.__dict__[fldname] = scale(PV)

        
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
    """
    def create_ijvecs(self):
        if not hasattr(self,'ns'):
            self.init_nasa()
        ns = self.ns
        kd = cKDTree(list(np.vstack((np.ravel(ns.lon),np.ravel(ns.lat))).T))
        dist,ij = kd.query(list(np.vstack((np.ravel(self.llon),
                                           np.ravel(self.llat))).T))
        l3j,l3i = np.meshgrid(np.arange(ns.j2-ns.j1),
                              np.arange(ns.i2-ns.i1))
        l3ij = np.vstack((np.ravel(l3i),np.ravel(l3j))).T
        gomj,gomi = np.meshgrid(np.arange(j2-j1),np.arange(i2-i1))
        #gomij = np.array(zip(np.ravel(gomi),np.ravel(gomj)))
        gomij = np.vstack((np.ravel(gomi),np.ravel(gomj))).T
        def gla(kk):
            mat = self.imat*0
            mat[self.ijvec[:,0],self.ijvec[:,1]] = l3ij[ij,:][:,kk]
            #for i,j in gomij:
            #    mat[i,j] = l3ij[ij,:][n,kk]
            #    n+=1
            return mat
        self.si = gla(1)
        self.satj = np.ravel(gla(0))
        self.sati = np.ravel(gla(1))
        
    def create_satijvecs(self):
        i1=self.i1; i2=self.i2;j1=self.j1; j2=self.j2
        if not hasattr(self,'ns'):
            self.init_nasa()
        ns = self.ns
        kd = cKDTree(list(np.vstack((np.ravel(self.llon),
                                     np.ravel(self.llat))).T))
        dist,ij = kd.query(list(np.vstack((np.ravel(ns.lon),
                                           np.ravel(ns.lat))).T))
        gomj,gomi = np.meshgrid(np.arange(ns.j2-ns.j1),
                                np.arange(ns.i2-ns.i1))
        l3j,l3i = np.meshgrid(np.arange(j2-j1),np.arange(i2-i1))
        l3ij  = np.vstack((np.ravel(l3i),np.ravel(l3j))).T
        gomij = np.vstack((np.ravel(gomi),np.ravel(gomj))).T
        def gla(kk):
            mat = gomj[:]*0
            mat[gomij[:,0],gomij[:,1]] = l3ij[ij,:][:,kk]
            return mat
        self.gcmj = np.ravel(gomi)
        self.gcmi = np.ravel(gomj)
        self.si = gla(1)
        self.satj = np.ravel(gla(0))
        self.sati = np.ravel(gla(1))
    """
        
    def add_landmask(self):
        """
        from scipy.stats import nanmean
        nc = pycdf.CDF(self.datadir + '/oscar_vel2005.nc')
        u = nc.var('u')[:,:,:,:1081]
        msk = np.squeeze(nanmean(u,axis=0))
        msk = msk*0
        msk[np.isnan(msk)]=1
        self.landmask = gmtgrid.convert(msk,self.gr)
        """
