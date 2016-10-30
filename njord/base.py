from __future__ import print_function

import os,sys
import datetime
import warnings
import ftplib
import urlparse
from datetime import datetime as dtm

import numpy as np
import pylab as pl
from scipy.spatial import cKDTree
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.animation as animation


from njord.utils import lldist
import requests
import gmtgrid

import config

try:
    import projmap
    USE_BASEMAP = True
except:
    USE_BASEMAP = False
try:
    import figpref
    USE_FIGPREF = True
except:
    USE_FIGPREF = False
try:
    import pyresample as pr
    USE_PYRESAMPLE = True
except:
    USE_PYRESAMPLE = False

    
dtm = datetime.datetime

class Grid(object):
    """Base class of njord for lat-lon gridded 2D or 3D data """
    def __init__(self, **kwargs):
        """Initalize an instance based on a config file"""
        self.projname = kwargs.get("projname", "%s.%s" %
                                   (self.__module__.split(".")[-1],
                                    type(self).__name__))
        preset_dict = config.load('njord', self.projname, kwargs)
        for key,val in preset_dict.items():
            setattr(self, key, val)
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)
        self.setup_grid()
        self._set_maxmin_ij()
                
    def _llboundaries_to_ij(self):
        """Calculate i1,i2,j1,j2 from lat and lon"""
        if hasattr(self, 'ijarea'):
            self.i1,self.i2,self.j1,self.j2 = getattr(self, 'ijarea')
        for att,val in zip(['i1', 'i2', 'j1', 'j2'], [0,self.imt,0,self.jmt]): 
            if not hasattr(self,att):
                self.__dict__[att] = val
        if hasattr(self, 'latlon'):
            self.lon1,self.lon2,self.lat1,self.lat2 = getattr(self, 'latlon')
        
        lat = self.llat if hasattr(self, 'llat') else self.latvec[:, np.newaxis]
        lon = self.llon if hasattr(self, 'llon') else self.lonvec[np.newaxis, :]
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

    @property
    def shape(self):
        return self.llat.shape
            
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
        if "date" in kwargs:
            self.jd = pl.datestr2num(kwargs['date'])
            self.jd = int(self.jd) if self.jd == int(self.jd) else self.jd
        elif ('yd' in kwargs) & ('yr' in kwargs):
            if self.yd < 1:
                self.yr = self.yr -1
                ydmax = (pl.date2num(dtm(self.yr, 12, 31)) -
                         pl.date2num(dtm(self.yr,  1,  1))) + 1    
                self.yd = ydmax + self.yd     
            self.jd = self.yd + pl.date2num(dtm(self.yr,1,1)) - 1
        elif  ('yr' in kwargs) & ('mn' in kwargs):
            if not "dy" in kwargs:
                kwargs["dy"] = 1
                setattr(self, "dy", 1)
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
        if self.jd < self.fulltvec.min():
            raise ValueError, "Date before first available model date"
        if self.jd > self.fulltvec.max():
            raise ValueError, "Date after last available model date"
        self._jd_to_dtm()


    def dx_approx(self):
        """Caclulate dx from llon and llat"""
        if not hasattr(self, "_dx_approx"):
            dx = self.llon * np.nan
            for i in xrange(self.jmt):
                latvec1 =  (self.llat[i,0:-2] + self.llat[i,1:-1]) / 2
                latvec2 =  (self.llat[i,2:]   + self.llat[i,1:-1]) / 2
                lonvec1 =  (self.llon[i,0:-2] + self.llon[i,1:-1]) / 2
                lonvec2 =  (self.llon[i,2:]   + self.llon[i,1:-1]) / 2
                dx[i,1:-1] = lldist.ll2dist2vec(lonvec1, latvec1,
                                                lonvec2, latvec2)
            dx[:, 0] = 2 * dx[:,1]  - dx[:,2]
            dx[:,-1] = 2 * dx[:,-2] - dx[:,-3]
            self._dx_approx = dx
        return self._dx_approx

    def dy_approx(self):
        """Caclulate dy from llon and llat"""
        if not hasattr(self, "_dx_approx"):
            dy = self.llat * np.nan
            for j in xrange(self.imt):
                latvec1 =  (self.llat[0:-2,j] + self.llat[1:-1,j]) / 2
                latvec2 =  (self.llat[2:,j]   + self.llat[1:-1,j]) / 2
                lonvec1 =  (self.llon[0:-2,j] + self.llon[1:-1,j]) / 2
                lonvec2 =  (self.llon[2:,j]   + self.llon[1:-1,j]) / 2
                dy[1:-1,j] = lldist.ll2dist2vec(lonvec1, latvec1,
                                                lonvec2, latvec2)
            dy[ 0,:] = 2 * dy[ 1,:] - dy[ 2,:]
            dy[-1,:] = 2 * dy[-2,:] - dy[-3,:]
            self._dy_approx = dy
        return self._dy_approx


    @property
    def datestr(self):
        jd = getattr(self, "jd", self.defaultjd)
        if type(jd) is int:
            return pl.num2date(jd).strftime("%Y-%m-%d")
        else:
            return pl.num2date(jd).strftime("%Y-%m-%d %H:%M")
    
        
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
            print(string)

                
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
            print(ipos, jpos)
            print(field[jpos, ipos], np.mean(field[jpos, ipos]))
            print(weights[dpos,:],  weights[dpos,:].sum())
            print(field[jpos, ipos] * weights[dpos,:])
            print((field[jpos, ipos] * weights[dpos,:]).sum())
            print(fldvec[dpos])
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
                print(intjdvec.max() - jd)
                mask = intjdvec == jd
                try:
                    fld = self.get_field(fldname, jd=jd)
                except IOError:
                    print("No file found")
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
    """
    def reproject(self, nj_obj, field):
        #Reproject a field of another njord inst. to the current grid
        if not hasattr(self,'nj_ivec'):
            self.add_njijvec(nj_obj)
        field = getattr(nj_obj, field) if type(field) is str else field
        
        if hasattr(nj_obj, 'tvec') and (len(nj_obj.tvec) == field.shape[0]):
            newfield = np.zeros(nj_obj.tvec.shape + self.llat.shape)
            for tpos in range(len(nj_obj.tvec)):
                newfield[tpos,:,:] = self.reproject(nj_obj, field[tpos,...])
            return newfield

        di = self.i2 - self.i1
        dj = self.j2 - self.j1
        xy = np.vstack((self.nj_jvec, self.nj_ivec))
        if type(field) == str:
            weights = np.ravel(nj_obj.__dict__[field])[self.nj_mask]
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
    """




    def reproject(self, nj_inst, field, roi=None):
        """Reproject a field of another njord inst. to the current grid"""
        field = getattr(nj_inst, field) if type(field) is str else field

        """
        if hasattr(nj_obj, 'tvec') and (len(nj_obj.tvec) == field.shape[0]):
            newfield = np.zeros(nj_obj.tvec.shape + self.llat.shape)
            for tpos in range(len(nj_obj.tvec)):
                newfield[tpos,:,:] = self.reproject(nj_obj, field[tpos,...])
            return newfield
        """
        
        if not hasattr(self, "_prdef"):
            self._prdef = pr.geometry.GridDefinition(lons=self.llon,
                                                     lats=self.llat)
        if not hasattr(nj_inst, "_prdef"):
            nj_inst._prdef = pr.geometry.GridDefinition(lons=nj_inst.llon,
                                                        lats=nj_inst.llat)
        if roi is None:
            roi = nj_inst.dy_approx().max() * 2

        return pr.kd_tree.resample_nearest(nj_inst._prdef, field,
                                           self._prdef, roi)
    
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
        return np.nanmedian(np.nanmedian(field.reshape(sh),axis=-1), axis=1)

    def get_field(self, field,  **kwargs):
        self.load(field, **kwargs)
        return self.__dict__[field]

    def retrive_file(self, url, local_filename=None, params=None):
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
            print("downloading\n %s \nto\n %s" % (url, local_filename))
            r = requests.get(url, params=params, stream=True,timeout=2)
            if r.ok:
                if local_filename is None:
                    return r.text
                else:
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

    @property
    def fulltvec(self):
        minjd = getattr(self, "minjd", 0)
        maxjd = getattr(self, "maxjd", int(pl.date2num(dtm.now())))
        return np.arange(minjd, maxjd+1)

    def get_tvec(self, jd1, jd2):
        jd1 = pl.datestr2num(jd1) if type(jd1) is str else jd1
        jd2 = pl.datestr2num(jd2) if type(jd2) is str else jd2
        tvec = self.fulltvec
        if jd1 < tvec.min():
            raise ValueError, "jd1 too small"
        if jd2 > tvec.max():
            raise ValueError, "jd2 too large"
        return tvec[(tvec >= jd1) & (tvec <= jd2)]
        
    def timeseries(self, fieldname, jd1, jd2, mask=None):
        """Create a timeseries of fields using mask to select data"""
        mask = mask if mask is not None else self.llat == self.llat
        self.tvec = self.get_tvec(jd1, jd2)
        field = np.zeros((len(self.tvec),) + self.llat.shape, dtype=np.float32)
        for n,jd in enumerate(self.tvec):
            print(pl.num2date(jd), len(self.tvec) - n)
            try:
                field[n,:,:] = self.get_field(fieldname, jd=jd).astype(np.float32)
            except (KeyError, IOError):
                field[n,:,:] = np.nan
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
            print(pl.num2date(jd))
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
    def mp(self, map_region=None):
        """Return handle to projmap instance"""
        if not hasattr(self, '_map_instance'):
            if USE_BASEMAP:
                self._map_instance = projmap.Projmap(self.map_region)
                self.mpxll,self.mpyll = self.mp(self.llon,self.llat)
            else:
                raise ImportError("Module 'projmap' not available")
        return self._map_instance

    def change_map_region(self, new_map_region):
            self.map_region = new_map_region
            if hasattr(self,'_map_instance'):
                    del self._map_instance


    def pcolor(self, fld, title=None, colorbar=False, **kwargs):
        """Make a pcolor-plot of field"""
        if USE_FIGPREF: figpref.current()
        if not hasattr(self, 'mp'): self.add_mp()
        miv = np.ma.masked_invalid
        if type(fld) == str:
            field = np.squeeze(getattr(self, fld))
        else:
            field = np.squeeze(fld)
        x,y = self.mp(self.llon, self.llat)
        self.mp.pcolormesh(x,y,miv(field), **kwargs)
        self.mp.nice()
        if title is not None:
            pl.title(title)        
        if colorbar:
            if type(colorbar) is dict:
                self.mp.colorbar(**colorbar)
            else:
                self.mp.colorbar()


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
        cs = self.mp.contour(x,y,miv(field), *args, **kwargs)
        self.mp.nice()
        return cs
        
    #    def movie(self,fld, jdvec, k=0, c1=0,c2=10):
    #        """Create a movie of a field """
    #        import anim
    #        mv = anim.Movie()
    #        for n,jd in enumerate(jdvec):
    #            self.load(fld,jd=jd)
    #            pl.clf()
    #            self.pcolor(self.__dict__[fld][k,:,:])
    #            pl.clim(c1, c2)
    #            pl.colorbar(pad=0,aspect=40)
    #            pl.title('%s %s' % (fld, pl.num2date(jd)
    #                                .strftime('%Y-%m-%d %H:%M')))
    #            mv.image()
    #        mv.video()



    def movie(self, fldname, jd1=None, jd2=None, jdvec=None, fps=10, **kwargs):
        curr_backend = plt.get_backend()
        plt.switch_backend('Agg')
        FFMpegWriter = animation.writers['ffmpeg']
        metadata = dict(title='%s' % (self.projname),
                        artist=self.projname,
                        comment='https://github.com/brorfred/njord')
        writer = FFMpegWriter(fps=fps, metadata=metadata,
            extra_args=['-vcodec', 'libx264',"-pix_fmt", "yuv420p"])

        jdvec = self.get_tvec(jd1, jd2) if jdvec is None else jdvec
        fig = plt.figure()
        with writer.saving(fig, "%s.mp4" % self.projname, 200):
            for jd in jdvec:
                pl.clf()
                print(pl.num2date(jd).strftime("%Y-%m-%d %H:%M load "), end="")
                sys.stdout.flush()
                try:
                    fld= self.get_field(fldname, jd=jd)
                except:
                    print("not downloaded" % jd)
                    continue
                print("plot ", end="")
                sys.stdout.flush()
                self.pcolor(fld, **kwargs)
                pl.title(pl.num2date(jd).strftime("%Y-%m-%d %H:%M"))
                print("write")
                writer.grab_frame()#bbox_inches="tight", pad_inches=0)
        plt.switch_backend(curr_backend)

    def check_gridtype(self):
        pass
