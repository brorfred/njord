
from __future__ import print_function

import os,sys
import datetime
import warnings
import ftplib
from urllib.parse import urlsplit
import shelve
import pathlib

import numpy as np
from scipy.spatial import cKDTree
import netCDF4
try:
    import pylab as pl
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import matplotlib.animation as animation
    HAS_MATPLOTLIB = False
except (ImportError, RuntimeError):
    import njord.utils.mpl_dates as pl
    HAS_MATPLOTLIB = False

import requests
import click

from njord.utils import lldist, yrday, time, gmtgrid
from njord import config

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
        if kwargs.get("mp") is not None:
            kwargs["map_region"] = kwargs["mp"]
        kwargs.pop("mp", None)
        if kwargs.get("map_region") is not None:
            mp = self._map_instance = projmap.Projmap(kwargs["map_region"])
            try:
                kwargs["lat1"] = kwargs.get("lat1", mp.llcrnrlat)
                kwargs["lat2"] = kwargs.get("lat2", mp.urcrnrlat)
                kwargs["lon1"] = kwargs.get("lon1", mp.llcrnrlon)
                kwargs["lon2"] = kwargs.get("lon2", mp.urcrnrlon)
            except:
                pass
             
        preset_dict = config.load('njord', self.projname, kwargs)
        for key,val in preset_dict.items():
            try:
                setattr(self, key, val)
            except AttributeError as err:
                print(err)
                raise AttributeError(key)
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)

        self.missing_fields = shelve.open(
            os.path.join(self.datadir, "missing_files.dbm"), writeback=True)

        self.setup_grid()
        if hasattr(self, "imt"):
            self._fullimt = self.imt
        if hasattr(self, "jmt"):
            self._fulljmt = self.jmt
        self._set_maxmin_ij()
        if (self._fullimt != self.imt) or (self._fulljmt != self.jmt):
            self._partial_grid = True
      
    def _llboundaries_to_ij(self):
        """Calculate i1,i2,j1,j2 from lat and lon"""
        if hasattr(self, 'ijarea'):
            self.i1,self.i2,self.j1,self.j2 = getattr(self, 'ijarea')
        for key,val in zip(['i1', 'i2', 'j1', 'j2'], [0,self.imt,0,self.jmt]): 
            if not hasattr(self, key):
                setattr(self, key, val)
        if hasattr(self, 'latlon'):
            self.lon1,self.lon2,self.lat1,self.lat2 = getattr(self, 'latlon')
        
        lat = self.llat if hasattr(self, 'llat') else self.latvec[:, np.newaxis]
        lon = self.llon if hasattr(self, 'llon') else self.lonvec[np.newaxis, :]
        if hasattr(self,'lat1'):
            if lat.min() >= self.lat1:
                self.lat1 = lat.min()
            elif lat[-1,0] < lat[0,0]:
                self.j2 = int(np.nonzero(lat>self.lat1)[0].max() + 1)
            else:
                self.j1 = int(np.nonzero(lat<self.lat1)[0].max())
        if hasattr(self,'lat2'):
            if lat.max() <= self.lat2:
                self.lat2 = lat.max()
            elif lat[-1,0] < lat[0,0]: 
                self.j1 = int(np.nonzero(lat<self.lat2)[0].min()-1)
            else:
                self.j2 = int(np.nonzero(lat>self.lat2)[0].min() + 1)
        if hasattr(self,'lon1'):
            if lon.min() >= self.lon1:
                self.lon1 = lon.min()
            else:
                self.i1 = int(np.nonzero(lon<=self.lon1)[-1].max()) 
        if hasattr(self,'lon2'):
            if lon.max() <= self.lon2:
                self.lon2 = lon.max()
            else:
                self.i2 = int(np.nonzero(lon>=self.lon2)[-1].min() + 1)
            
    @property
    def shape(self):
        """Return the shape of the field"""
        return self.llat.shape
            
    def _set_maxmin_ij(self):
        """Set a potential sub region of the grid if defined """
        self._llboundaries_to_ij()
        if hasattr(self, "latvec"):
            self.latvec = self.latvec[self.j1:self.j2]
        if hasattr(self, "lonvec"):
            self.lonvec = self.lonvec[self.i1:self.i2]
        if self.flipped_y:
            i2 = None if self.jmt == self.j2 else self.jmt-self.j2
            self.ijslice_flip = (slice(self.jmt-self.j1, i2, -1),
                                 slice(         self.i1,          self.i2,  1))
        self.ijslice = (slice(self.j1,self.j2), slice(self.i1,self.i2))
        for var in ['depth', 'dxt', 'dyt',  'dxu', 'dyu',  'dzt',
                    'area',  'vol',  'llon', 'llat']:
            if hasattr(self,var):
                setattr(self, var, getattr(self, var)[self.ijslice])
        self.imt = self.i2 - self.i1
        self.jmt = self.j2 - self.j1
        
    def _timeparams(self, **time_kwargs):
        """Calculate time parameters from given values"""
        tdict = time.expand_timeparams(self, **time_kwargs)
        for key in tdict:
            setattr(self, key, tdict[key])
        if self.jd < self.minjd:
            raise ValueError("Date before first available model date")
        if self.jd > self.maxjd:
            raise ValueError("Date after last available model date")

    def modeljd(self, jd):
        return self.fulljdvec[self.fulljdvec <= jd].max()

    def _model_timeparams(self, **time_kwargs):
        tdict = time.expand_timeparams(self, **time_kwargs)
        jd = self.modeljd(tdict["jd"])
        self._timeparams(jd=jd)
        
    def dx_approx(self):
        """Calculate dx from llon and llat"""
        if not hasattr(self, "_dx_approx"):
            dx = self.llon * np.nan
            for i in range(self.jmt):
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
        if not hasattr(self, "_dy_approx"):
            dy = self.llat * np.nan
            for j in range(self.imt):
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
    def datadir(self):
        return self._datadir

    @datadir.setter
    def datadir(self, path):
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)
        self._datadir = path
    
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
        """Generate a KD-tree objects for the current njord instance"""
        if coord == "ij":
            self.add_ij()
            xvec = self.kdivec = np.ravel(self.imat)
            yvec = self.kdjvec = np.ravel(self.jmat)
        else:
            xvec = self.kdlonvec = self.llon.flat
            yvec = self.kdlatvec = self.llat.flat
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
        
    def ll2ij(self, lonvec, latvec, mask=None, cutoff=None,
                  nei=1, all_nei=False, return_dist=False):
        """Reproject a lat-lon vector to i-j grid coordinates"""
        self.add_ij()
        if mask is not None:
            self.add_kd(mask)
        elif not hasattr(self, 'kd'):
            self.add_kd()
        dist,ij = self.kd.query(list(np.vstack((lonvec, latvec)).T), nei)
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

    def isncfile(self, filename):
        if not os.path.isfile(filename):
            return False
        try:
            netCDF4.Dataset(filename)
            return True
        except:
            os.unlink(filename)
            return False

    def retrive_file(self, url, local_filename=None,username=None,params=None):
        """Retrive file from remote server via http

        https://stackoverflow.com/questions/51684008
        """
        spliturl = urlsplit(url)
        self.vprint(spliturl)
        #spliturl = urllib.parse.urlsplit(url)
        if spliturl.scheme == "ftp":
            try:
                ftp = ftplib.FTP(spliturl.netloc) 
                self.vprint(ftp.login(self.ftpuser, self.ftpasswd))
                ftpdir = spliturl.path
                self.vprint("Change dir to '%s'" % ftpdir)
                self.vprint(ftp.cwd(ftpdir))
            except ftplib.error_perm as err:
                print (spliturl.netloc)
                print (os.path.split(spliturl.path)[0])
                raise IOError(err)

            #self.vprint(local_filename)
            lfn = os.path.basename(local_filename)
            if not lfn in ftp.nlst():
                print(ftp.nlst())
                raise ftplib.Error("'%s' is not the ftp server" % lfn)
            with open(local_filename, 'wb') as lfh:
                ftp.voidcmd('TYPE I')
                length = ftp.size(lfn)
                short_lfn = lfn if len(lfn)<18 else lfn[:9] + "..." + lfn[-9:]
                with click.progressbar(length=length, label=short_lfn) as bar:
                    def file_write(data):
                        lfh.write(data) 
                        bar.update(len(data))
                    try:
                        ftp.retrbinary("RETR %s" % lfn, file_write)
                    except ftplib.error_perm as err:
                        os.unlink(local_filename)
                        raise IOError(err)
            ftp.quit()
            return True
        else:
            print("downloading\n %s \nto\n %s" % (url, local_filename))
            try:
                r = requests.get(url, params=params, stream=True, timeout=2)
            except requests.ReadTimeout:
                warnings.warn("Connection to server timed out.")
                return False
            if r.ok:
                if local_filename is None:
                    return r.text
                else:
                    with open(local_filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=1024): 
                            if chunk:
                                f.write(chunk)
                                f.flush()
                    return True
            else:
                warnings.warn("Could not download file from server")
                return False

    def fatten_landmask(self):
        if not hasattr(self,'landmask'): self.add_landmask()
        i,j = np.nonzero(self.landmask==0)
        ii = i[(i>0) & (j>0) & (i<self.imt-1) & (j<self.jmt-1)]
        jj = j[(i>0) & (j>0) & (i<self.imt-1) & (j<self.jmt-1)]
        for i,j in zip([0,0,1,-1, -1,-1,1,1],[1,-1,0,0,-1,1,-1,1]):
            self.landmask[ii+i,jj+j]=0

    @property
    def fulljdvec(self):
        if self.fldsperday != 1:
            fpd = self.fldsperday
            return np.arange(self.minjd*fpd, self.maxjd*fpd + 1)/fpd
        else:
            dpf = self.daysperfld
        return np.arange(self.minjd, self.maxjd+dpf, dpf)

    @property
    def yrvec(self):
        return yrday.years(self.fulljdvec)

    @property
    def mnvec(self):
        return yrday.months(self.fulljdvec)
    
    def get_tvec(self, jd1, jd2):
        jd1 = pl.datestr2num(jd1) if type(jd1) is str else jd1
        jd2 = pl.datestr2num(jd2) if type(jd2) is str else jd2
        tvec = self.fulljdvec
        if jd1 < tvec.min():
            raise ValueError("jd1 too small")
        if jd2 > tvec.max():
            raise ValueError("jd2 too large")
        return tvec[(tvec >= jd1) & (tvec <= jd2)]
        
    def timeseries(self, fieldname, jd1, jd2, mask=None, **loadkwargs):
        """Create a timeseries of fields using mask to select data"""
        mask = mask if mask is not None else self.llat == self.llat
        self.tvec = self.get_tvec(jd1, jd2)
        field = np.zeros((len(self.tvec),) + self.llat.shape, dtype=np.float32)
        for n,jd in enumerate(self.tvec):
            print(pl.num2date(jd), len(self.tvec) - n)
            try:
                field[n,:,:] = self.get_field(
                    fieldname, jd=jd, **loadkwargs).astype(np.float32)
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
        im = self.mp.pcolormesh(x,y,miv(field), **kwargs)
        self.mp.nice()
        if title is not None:
            pl.title(title)        
        if colorbar:
            if type(colorbar) is dict:
                self.mp.colorbar(**colorbar)
            else:
                self.mp.colorbar()
        return im

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
