import os
import ConfigParser
import json

import numpy as np
from scipy.spatial import cKDTree

import projmap

class Grid(object):
    """Base class of njord for lat-lon gridded 2D or 3D data """
    def __init__(self, **kwargs):
        """Initalize an instance based on a config file"""
        self.projname = "%s.%s" % (self.__module__.split(".")[-1],
                                   type(self).__name__)
        self.basedir =  os.path.dirname(os.path.abspath(__file__))
        self.inkwargs = kwargs
        self._load_presets('njord',kwargs)
        for key,val in kwargs.items():
            self.__dict__[key] = val
        self.setup_grid()
        self._set_maxmin_ij()

    def _load_presets(self, filepref, kwargs):
        """Read and parse the config file"""
        cfg = ConfigParser.ConfigParser()
        files = ["%s/%s.cfg" % (os.curdir, filepref),
                 "%s/.%s.cfg" % (os.path.expanduser("~"), filepref),
                 "%s/%s.cfg" % (self.basedir, filepref)]
        for fnm in files:
            cfg.read(fnm)
            if self.projname in cfg.sections():
                self.config_file = fnm
                break
        else:
            raise NameError('Project not included in config files')

        def splitkey(key, val):
            if key in self.inkwargs.keys():
                self.__dict__[key] = self.inkwargs[key]
                del self.inkwargs[key]
            else:
                self.__dict__[key] = val
            
        for key,val in cfg.items(self.projname):
            try:
                splitkey(key, json.loads(val))
            except ValueError:
                splitkey(key, val)


    def _set_maxmin_ij(self):
        """Set a potential sub region of the grid if defined """
        if 'ijarea' in self.inkwargs.keys():
            self.i1,self.i2,self.j1,self.j2 = self.inkwargs['ijarea']    
        for att,val in zip(['i1', 'i2', 'j1', 'j2'], [0,self.imt,0,self.jmt]): 
            if not hasattr(self,att):
                self.__dict__[att] = val

        if 'latlon' in self.inkwargs.keys():
            self.lon1,self.lon2,self.lat1,self.lat2 = self.inkwargs['latlon']
        if hasattr(self,'lat1'):
            self.j1 = np.nonzero(self.llat<self.lat1)[0].max() - 1
        if hasattr(self,'lat2'):
            self.j2 = np.nonzero(self.llat>self.lat2)[0].min() + 1
        if hasattr(self,'lon1'):
            self.i1 = np.nonzero(self.llon<self.lon1)[1].max() - 1 
        if hasattr(self,'lon2'):
            self.i2 = np.nonzero(self.llon>self.lon2)[1].min() + 1

        self.llat = self.llat[self.j1:self.j2, self.i1:self.i2]
        self.llon = self.llon[self.j1:self.j2, self.i1:self.i2]

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
        """Create a timeseries of fields using mask to select data"""
        if len(mask)==0: mask = (self.lat !=-800)
        field = self.llat * 0 
        for n,jd in enumerate(np.arange(jd1,jd2+1)):
            self.load(fieldname, jd=jd)
            field[n,:,:] = self.__dict__[fieldname]
            field[n,~mask] = np.nan
        self.__dict__[fieldname + 't'] = field
        self.tvec = np.arange(jd1, jd2+1)

    def histmoller(self, fieldname, jd1, jd2, y1, y2,
                   mask=[], bins=100, logy=True):
        """Create a histmoller diagram (histogram on y and time on x)"""
        if len(mask)==0: mask = (self.llat !=-800)
        if logy:
            vlist = np.exp(np.linspace(np.log(y1), np.log(y2), bins+1))
        else:
            vlist = np.linspace(y1, y2, bins+1)
        hsmat = np.zeros((jd2-jd1+1,bins), dtype=np.int)
        tvec = np.arange(jd1,jd2+1)
        for n_t,jd in enumerate(tvec):
            print pl.num2date(jd)
            self.load(fieldname, jd=jd)
            field = self.__dict__[fieldname]
            field[~mask] = np.nan
            hsmat[n_t,:],_ = np.histogram(field[~np.isnan(field)], vlist)
        class hiS: pass
        hiS.tvec = tvec
        hiS.vlist = vlist
        hiS.hist = hsmat
        hiS.norm = (hsmat.astype(float) / hsmat.max(axis=0))
        self.__dict__[fieldname + "_hiS"] = hiS

    def dump(self, filename=None):
        """Dump all attributes in instance to a file"""
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
        self.mp = projmap.Projmap(self.map_region)
        
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
