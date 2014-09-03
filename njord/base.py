import os
import ConfigParser
import json
import datetime
import warnings
import ftplib
import urlparse

import numpy as np
import pylab as pl
from scipy.spatial import cKDTree
from scipy.stats import nanmedian

import requests
import gmtgrid
import projmap
try:
    import figpref
    USE_FIGPREF = True
except:
    USE_FIGPREF = False

dtm = datetime.datetime

class Grid(object):
    """Base class of njord for lat-lon gridded 2D or 3D data """
    def __init__(self, **kwargs):
        """Initalize an instance based on a config file"""
        if not hasattr(self, 'projname'):
            self.projname = "%s.%s" % (self.__module__.split(".")[-1],
                                       type(self).__name__)
        self.basedir =  os.path.dirname(os.path.abspath(__file__))
        self.inkwargs = kwargs
        self._load_presets('njord',kwargs)
        for key,val in kwargs.items():
            setattr(self, key, val)
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)
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

    def _llboundaries_to_ij(self):
        """Calculate i1,i2,j1,j2 from lat and lon"""
        if 'ijarea' in self.inkwargs.keys():
            self.i1,self.i2,self.j1,self.j2 = self.inkwargs['ijarea']    
        for att,val in zip(['i1', 'i2', 'j1', 'j2'], [0,self.imt,0,self.jmt]): 
            if not hasattr(self,att):
                self.__dict__[att] = val
        if 'latlon' in self.inkwargs.keys():
            self.lon1,self.lon2,self.lat1,self.lat2 = self.inkwargs['latlon']
        
        lat = self.llat if hasattr(self, 'llat') else self.latvec[:, np.newaxis]
        lon = self.llon if hasattr(self, 'llon') else self.lonvec[:, np.newaxis]
        if hasattr(self,'lat1'):
            if lat[-1,0] < lat[0,0]:
                self.j2 = int(np.nonzero(lat>self.lat1)[0].max() + 1)
            else:
                self.j1 = int(np.nonzero(lat<self.lat1)[0].max())
        if hasattr(self,'lat2'):
            if lat[-1,0] < lat[0,0]: 
                self.j1 = int(np.nonzero(lat<self.lat2)[0].min()-1)
            else:
                self.j2 = int(np.nonzero(lat>self.lat2)[0].min() + 1)
        if hasattr(self,'lon1'):
            self.i1 = int(np.nonzero(lon<=self.lon1)[-1].max()) 
        if hasattr(self,'lon2'):
            self.i2 = int(np.nonzero(lon>=self.lon2)[-1].min() + 1) 

    def _set_maxmin_ij(self):
        """Set a potential sub region of the grid if defined """
        self._llboundaries_to_ij()
        if hasattr(self, "latvec"):
            self.latvec = self.latvec[self.j1:self.j2]
        if hasattr(self, "lonvec"):
            self.lonvec = self.lonvec[self.i1:self.i2]                
        for var in ['depth', 'dxt', 'dyt',  'dxu', 'dyu',  'dzt',
                    'area',  'vol',  'llon', 'llat']:
            if hasattr(self,var):
                self.__dict__[var] = self.__dict__[var][self.j1:self.j2,
                                                        self.i1:self.i2]

    def _timeparams(self, **kwargs):
        """Calculate time parameters from given values"""
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
        if ('yd' in kwargs) & ('yr' in kwargs):
            if self.yd < 1:
                self.yr = self.yr -1
                ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                         pl.date2num(dtm(self.yr,  1,  1))) + 1    
                self.yd = ydmax + self.yd     
            self.jd = self.yd + pl.date2num(dtm(self.yr,1,1)) - 1
        elif  ('yr' in kwargs) & ('mn' in kwargs) & ('dy' in kwargs):
            self.jd = pl.date2num(dtm(self.yr,self.mn,self.dy))
        elif not 'jd' in kwargs:
            if hasattr(self, 'defaultjd'):
                self.jd = self.defaultjd
            else:
                raise KeyError, "Time parameter missing"
        if hasattr(self,'hourlist'):
            dd = self.jd-int(self.jd)
            ddlist = np.array(self.hourlist).astype(float)/24
            ddpos = np.argmin(np.abs(ddlist-dd))
            self.jd = int(self.jd) + ddlist[ddpos]
        self._jd_to_dtm()

    def _jd_to_dtm(self):
        dtobj = pl.num2date(self.jd)
        njattrlist = ['yr',  'mn',   'dy', 'hr',  'min',    'sec']
        dtattrlist = ['year','month','day','hour','minute', 'second']
        for njattr,dtattr in zip(njattrlist, dtattrlist):
            setattr(self, njattr, getattr(dtobj, dtattr))
        self.yd = self.jd - pl.date2num(dtm(self.yr,1,1)) + 1

    def _calcload(self, fldname, **kwargs):
        if fldname == "uv":
            self.load('u', **kwargs)
            self.load('v', **kwargs)
            km,jm,im = np.minimum(self.u.shape, self.v.shape)
            self.uv = np.sqrt(self.u[:km, :jm, :im]**2 + 
                              self.v[:km, :jm, :im]**2)/2
            return True
        else:
            return False


    def vprint(self, string, log_level=None):
        if getattr(self, 'verbose', False) == True:
            print string

                
    def add_ij(self):
        self.imat,self.jmat = np.meshgrid(np.arange(self.i2-self.i1),
                                          np.arange(self.j2-self.j1))
        self.kdijvec = np.vstack((np.ravel(self.imat),
                                  np.ravel(self.jmat))).T

        self._ijvec = np.arange(np.prod(self.imat.shape))

    @property
    def ijvec(self):
        if not hasattr(self, '_ijvec'):
            self.add_ij()
        return self._ijvec
            
    def add_kd(self, mask=None, coord="lonlat"):
        """Generate a KD-tree objects and for the current njord instance"""
        if coord == "ij":
            self.add_ij()
            xvec = self.kdivec = np.ravel(self.imat)
            yvec = self.kdjvec = np.ravel(self.jmat)
        else:
            xvec = self.kdlonvec = np.ravel(self.llon)
            yvec = self.kdlatvec = np.ravel(self.llat)
        if not mask is None: 
            xvec = xvec[~np.isnan(np.ravel(mask))]
            yvec = yvec[~np.isnan(np.ravel(mask))]
            self.add_ij()
            self.kdijvec = self.kdijvec[~np.isnan(np.ravel(mask))]
        self.kd = cKDTree(list(np.vstack((xvec, yvec)).T))

    def ijinterp(self,ivec,jvec, field, mask=None, nei=3, dpos=None):
        self.add_ij()
        if mask is not None:
            self.add_kd(mask,coord="ij")
        elif not hasattr(self, 'kdivec'):
            self.add_kd(coord="ij")
        dist,ij = self.kd.query(list(np.vstack((ivec,jvec)).T), nei)
        sumvec = np.nansum(field[self.kdijvec[ij][:,:,1],
                              self.kdijvec[ij][:,:,0]]*0+1,axis=1)
        weights = (1-dist/sumvec[:,np.newaxis])
        weights = weights/weights.sum(axis=1)[:,np.newaxis]
        fldvec = weights[:,0] * 0
        for n in np.arange(nei):
            fldvec = fldvec + field[self.kdijvec[ij[:,n]][:,1],
                                    self.kdijvec[ij[:,n]][:,0]] * weights[:,n]
        if dpos is not None:
            ipos = self.kdijvec[ij[dpos,:], 0]
            jpos = self.kdijvec[ij[dpos,:], 1]
            print ipos, jpos
            print field[jpos, ipos], np.mean(field[jpos, ipos])
            print weights[dpos,:],  weights[dpos,:].sum()
            print field[jpos, ipos] * weights[dpos,:]
            print (field[jpos, ipos] * weights[dpos,:]).sum()
            print fldvec[dpos]
        return fldvec
        
    def ll2ij(self, lon, lat, mask=None, cutoff=None, nei=1, all_nei=False,
              return_dist=False):
        """Reproject a lat-lon vector to i-j grid coordinates"""
        self.add_ij()
        if mask is not None:
            self.add_kd(mask)
        elif not hasattr(self, 'kd'):
            self.add_kd()
        dist,ij = self.kd.query(list(np.vstack((lon,lat)).T), nei)
        #if cutoff is not None:
        #    ij[dist>cutoff] = 0
        if nei == 1 :
            ivec = self.kdijvec[ij[:] - 1][:, 0]
            jvec = self.kdijvec[ij[:] - 1][:, 1]
            if cutoff is not None:
                ivec[dist>cutoff] = -999
                jvec[dist>cutoff] = -999
        elif all_nei == False:
            ivec = np.squeeze(self.kdijvec[ij[:,:]-1])[:, nei-1, 0]
            jvec = np.squeeze(self.kdijvec[ij[:,:]-1])[:, nei-1, 1]
            dist = np.squeeze(dist[:, nei-1])
        else:
            ivec = np.squeeze(self.kdijvec[ij[:,:]-1])[:, :, 0]
            jvec = np.squeeze(self.kdijvec[ij[:,:]-1])[:, :, 1]
        if return_dist == True:
            return ivec,jvec,dist
        else:
            return ivec,jvec

    def fld2vec(self, fldname, lonvec, latvec, jdvec, maskvec=None, djd=1,
                daysback=1, nei=1,all_nei=True):
        """Get data from the lat-lon-jd positions """
        if maskvec is not None:
            lonvec = lonvec[maskvec]
            latvec = latvec[maskvec]
            jdvec  = jdvec[maskvec]
        ivec,jvec = self.ll2ij(lonvec, latvec, nei=nei, all_nei=all_nei) 
        intjdvec  = ((jdvec / djd).astype(np.int)) * djd
        fldvec    = np.zeros((daysback, len(latvec))) * np.nan
        for days in np.arange(daysback):
            for n,jd in enumerate(np.unique(intjdvec)):
                print intjdvec.max() - jd
                mask = intjdvec == jd
                try:
                    fld = self.get_field(fldname, jd=jd)
                except IOError:
                    print "No file found"
                    continue
                #fldvec[days, mask] = self.ijinterp(ivec[mask],jvec[mask], fld)
                if ivec.ndim == 2:
                    fldvec[days, mask] = np.nanmean(fld[jvec, ivec], axis=1)[mask]
                else:
                    fldvec[days, mask] = fld[jvec[mask], ivec[mask]]
        return np.squeeze(fldvec)
            
    def add_ij2ij(self, njord_obj):
        sci,scj = njord_obj.ll2ij(np.ravel(self.llon),np.ravel(self.llat))
        self.sci = sci.reshape(self.i2,self.j2)
        self.scj = scj.reshape(self.i2,self.j2)

    def add_njijvec(self, njord_obj):
        lonvec = np.ravel(njord_obj.llon)
        latvec = np.ravel(njord_obj.llat)
        lonmsk = (lonvec>=self.llon.min()) & (lonvec<=self.llon.max())
        latmsk = (latvec>=self.llat.min()) & (latvec<=self.llat.max()) 
        self.nj_mask = lonmsk & latmsk
        self.nj_ivec,self.nj_jvec = self.ll2ij(lonvec[self.nj_mask],
                                               latvec[self.nj_mask])

    def reproject(self, njord_obj, field):
        """Reproject a field of another njord inst. to the current grid"""
        if not hasattr(self,'nj_ivec'): self.add_njijvec(njord_obj)
        di = self.i2 - self.i1
        dj = self.j2 - self.j1
        xy = np.vstack((self.nj_jvec, self.nj_ivec))
        if type(field) == str:
            weights = np.ravel(njord_obj.__dict__[field])[self.nj_mask]
        else:
            weights = np.ravel(field)[self.nj_mask]
        mask = ~np.isnan(weights) 
        flat_coord = np.ravel_multi_index(xy[:,mask],(dj, di))
        sums = np.bincount(flat_coord, weights[mask])
        cnts = np.bincount(flat_coord)
        fld = np.zeros((dj, di)) * np.nan
        fld.flat[:len(sums)] = sums.astype(np.float)/cnts
        try:
            self.add_landmask()
            fld[self.landmask] = np.nan
        except:
            print "Couldn't load landmask for %s" % self.projname
        return fld

    def decheck(self,fieldname=None):
        """Remove checkerboarding by smearing the field"""
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

    def rebin(self, field, shape):
        """Rebin field to a coarser matrix"""
        sh = shape[0],field.shape[0]//shape[0],shape[1],field.shape[1]//shape[1]
        return nanmedian(nanmedian(field.reshape(sh),axis=-1), axis=1)

    def get_field(self, field,  **kwargs):
        self.load(field, **kwargs)
        return self.__dict__[field]

    def retrive_file(self, url, local_filename):
        """Retrive file from remote server via http"""
        spliturl = urlparse.urlsplit(url)
        if spliturl.scheme == "ftp":
            ftp = ftplib.FTP(spliturl.netloc) 
            ftp.login("anonymous", "njord@bror.us") 
            ftp.cwd(spliturl.path)
            ftp.retrbinary("RETR %s" % os.path.basename(local_filename),
                           open(local_filename, 'wb').write)
            ftp.quit()
        else:
            r = requests.get(url, stream=True)
            if r.ok:
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024): 
                        if chunk: # filter out keep-alive new chunks
                            f.write(chunk)
                            f.flush()
                return True
            else:
                warnings.warn("Could not download file from server")
                return False

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

    def timeseries(self, fieldname, jd1, jd2, mask=None):
        """Create a timeseries of fields using mask to select data"""
        if len(mask)==0: mask = (self.lat !=-800)
        mask = mask if mask is not None else self.llat == self.llat
        for n,jd in enumerate(np.arange(jd1,jd2+1)):
            self.load(fieldname, jd=jd)
            field[n,:,:] = self.__dict__[fieldname]
            field[n,~mask] = np.nan
        self.__dict__[fieldname + 't'] = field
        self.tvec = np.arange(jd1, jd2+1)
        field = np.zeros((len(self.tvec),) + self.llat.shape) 
        for n,jd in enumerate(self.tvec):
            print pl.num2date(jd), pl.num2date(jd2)
            field[n,:,:] = self.get_field(fieldname, jd=jd)
            field[n, ~mask] = np.nan
        setattr(self, fieldname + 't', field)

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


    def add_mp(self, map_region=None):
        """Old method, don't use!"""
        return

    @property
    def mp(self, **kwargs):
        """Return a projmap instance as defined by self.map_region"""
        if not 'mp' in self.__dict__.keys():
            self.__dict__['mp'] = projmap.Projmap(self.map_region, **kwargs)
        if self.__dict__['mp'].region != self.map_region:
            self.__dict__['mp'] = projmap.Projmap(self.map_region)
        return self.__dict__['mp']

    def pcolor(self,fld, **kwargs):
        """Make a pcolor-plot of field"""
        if USE_FIGPREF: figpref.current()
        if not hasattr(self, 'mp'): self.add_mp()
        miv = np.ma.masked_invalid
        if type(fld) == str:
            field = self.__dict__[fld]
        else:
            field = fld
        x,y = self.mp(self.llon, self.llat)
        self.mp.pcolormesh(x,y,miv(field), **kwargs)
        self.mp.nice()

    def contour(self,fld, *args, **kwargs):
        """Make a pcolor-plot of field"""
        if USE_FIGPREF: figpref.current()
        if not hasattr(self, 'mp'): self.add_mp()
        miv = np.ma.masked_invalid
        if type(fld) == str:
            field = self.__dict__[fld]
        else:
            field = fld
        x,y = self.mp(self.llon, self.llat)
        self.mp.contour(x,y,miv(field), *args, **kwargs)
        self.mp.nice()

        
    def movie(self,fld, jdvec, k=0, c1=0,c2=10):
        """Create a movie of a field """
        import anim
        mv = anim.Movie()
        for n,jd in enumerate(jdvec):
            self.load(fld,jd=jd)
            pl.clf()
            self.pcolor(self.__dict__[fld][k,:,:])
            pl.clim(c1, c2)
            pl.colorbar(pad=0,aspect=40)
            pl.title('%s %s' % (fld, pl.num2date(jd)
                                .strftime('%Y-%m-%d %H:%M')))
            mv.image()
        mv.video()

    def check_gridtype(self):
        pass
