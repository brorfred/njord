import numpy as np
from scipy.spatial import cKDTree

import projmap

class Grid(object):
    """Base class of njord for lat-lon gridded 2D or 3D data """
    def __init__(self, **kwargs):
        
        self.class_name = type(self).__name__
        self.module_name = self.__module__
        
        


    def add_ij(self):
        self.jmat,self.imat = np.meshgrid(np.arange(self.j2-self.j1),
                                          np.arange(self.i2-self.i1))
        self.kdijvec = np.vstack((np.ravel(self.imat),
                                  np.ravel(self.jmat))).T

    def add_kd(self,mask=None):
        self.kdlatvec = np.ravel(self.llat)
        self.kdlonvec = np.ravel(self.llon)
        if not mask is None: 
            self.kdlatvec = self.kdlatvec[~np.isnan(np.ravel(mask))]
            self.kdlonvec = self.kdlonvec[~np.isnan(np.ravel(mask))]
            self.add_ij()
            self.kdijvec = self.kdijvec[~np.isnan(np.ravel(mask))]
        self.kd = cKDTree(list(np.vstack((self.kdlonvec,
                                          self.kdlatvec)).T))
                                          
    def ll2ij(self,lon,lat,nei=1):
        if not hasattr(self, 'kd'):
            self.add_kd()
            self.add_ij()
        dist,ij = self.kd.query(list(np.vstack((lon,lat)).T), nei)
        if nei == 1 :
            return self.kdijvec[ij[:] - 1][:, 0], self.kdijvec[ij[:]-1][:, 1]
        else:
            return (np.squeeze(self.kdijvec[ij[:,:]-1])[:, nei-1, 0],
                    np.squeeze(self.kdijvec[ij[:,:]-1])[:, nei-1, 1])

    def add_ij2ij(self,sc):
        sci,scj = sc.ll2ij(np.ravel(self.llon),np.ravel(self.llat))
        self.sci = sci.reshape(self.i2,self.j2)
        self.scj = scj.reshape(self.i2,self.j2)

    def decheck(self,fieldname=None):
        if not fieldname: 
            fieldname = self.last_loaded_feld
        fld = self.__dict__[fieldname]
        fldi = (fld[:,:-1,:] + fld[:,1:,:]) / 2
        fldj = (fld[:,:,:-1] + fld[:,:,1:]) / 2
        flint = (fld[:,-1,:] + fld[:,-2,:]) / 2
        flint = flint[:,1:] + flint[:,:-1]
        fljnt = (fld[:,:,-1] + fld[:,:,-2]) / 2
        fljnt = fljnt[:,1:] + fljnt[:,:-1]
        fld[:,:-1,:-1] = (fldj[:,:-1,:] + fldi[:,:,:-1])/2
        fld[:,-1,:-1] = flint
        fld[:,:-1,-1] = fljnt

    def get_field(self, field, jd):
        self.load(field,jd=jd)
        return self.__dict__[field]

    def get_landmask(self):
        if not hasattr(self,'landmask'):
            self.add_landmask()
        return self.landmask

    def fatten_landmask(self):
        if not hasattr(self,'landmask'): self.add_landmask()
        i,j = np.nonzero(self.landmask==0)
        ii = i[(i>0) & (j>0) & (i<self.imt-1) & (j<self.jmt-1)]
        jj = j[(i>0) & (j>0) & (i<self.imt-1) & (j<self.jmt-1)]
        for i,j in zip([0,0,1,-1, -1,-1,1,1],[1,-1,0,0,-1,1,-1,1]):
            self.landmask[ii+i,jj+j]=0


    def timeseries(self, fieldname, jd1, jd2, mask=[]):
        if len(mask)==0: mask = (self.lat !=-800)
        field = np.zeros( (jd2-jd1+1,self.i2-self.i1, self.j2-self.j1),
                          dtype=np.float32)
        for n,jd in enumerate(np.arange(jd1,jd2+1)):
            self.load(fieldname, jd=jd)
            field[n,:,:] = self.__dict__[fieldname]
            field[n,~mask] = np.nan
        self.__dict__[fieldname + 't'] = field
        self.tvec = np.arange(jd1, jd2+1)

    def llrect2ij(self,lon1,lon2,lat1,lat2):
        i1 = np.nonzero(self.llat[:,0]>lat2)[0].max()
        i2 = np.nonzero(self.llat[:,0]<lat1)[0].min()
        j1 = np.nonzero(self.llon[0,:]<lon1)[0].max()
        j2 = np.nonzero(self.llon[0,:]>lon2)[0].min()
        return int(i1),int(i2),int(j1),int(j2)

    def dump(self, filename=None):
        """Dump all atrributes in instance to a file"""
        dmp_dct = {}
        for v in self.__dict__.keys():
            dmp_dct[v] = self.__dict__[v]
        np.savez('arfields.npz',**dmp_dct)

    def pickup(self, filename=None):
        """Pickup a dumped attribute file"""
        dmpfile = np.load('arfields.npz')
        for v in self.__dict__.keys():
            self.__dict__[v] = dmpfile[v]

    def add_mp(self):
        self.mp = projmap.Projmap(self.region)
        

    def movie(self,fld='temp', k=39,jd1=730120,jd2=730120+365):
        """Create movie of a field """
        import anim
        miv = np.ma.masked_invalid
        x,y = self.mp(self.llon, self.llat)

        mv = anim.Movie()
        for n,jd in enumerate(np.arange(jd1,jd2)):
            self.load(fld,jd)
            pl.clf()
            self.mp.pcolormesh(x,y,miv(self.__dict__[fld][k,:,:]))
            self.mp.nice()
            pl.clim(0,2)
            pl.colorbar(pad=0,aspect=40)
            pl.title('%s %s' % (fld, pl.num2date(jd).strftime('%Y-%m-%d')))
            mv.image()
        mv.video()
