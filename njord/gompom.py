import glob
from datetime import datetime as dtm

import numpy as np
from scipy.io import netcdf_file

from pyhdf.SD import SD,SDC

import base

class Casco(base.Grid):

    def __init__(self, **kwargs):
        super(Casco, self).__init__()
        self.add_mp()
        self.missing = {'temp':-25184,'salt':-26128,'elev':16608,
                        'u':512,'v':512}

    def setup_grid(self):
        gc = netcdf_file(self.gridfile)
        self.llat = gc.variables['y'][:].copy()
        self.llon = gc.variables['x'][:].copy()
        self.depth = gc.variables['depth'][:].copy()
        self.depth[self.depth<0] = np.nan
        
    def load(self,field,**kwargs):
        """ Load casco fields for a given day"""
        self.last_loaded_feld = field
        self._timeparams(**kwargs)        
        t  = int(np.floor(np.modf(self.jd)[0]*8))
        filename = ps.path.join(self.datadir,
                             "casco.%04i%02i%02i.cdf" (self.yr,self.mn,self.dy)
        nc = netcdf_file(filename)        

        if field == 'uv':
            self.load('u', self.jd)
            self.load('v', self.jd)
            ang = nc.var('ang')[:] / 180 * np.pi
            uu = ( self.u * np.cos(ang) + self.v * np.sin(ang) ) 
            self.v = (-self.u * np.sin(ang) + self.v * np.cos(ang) )
            self.u = uu
        else:
            var = nc.variables[field]
            fld = var[t, self.i1:self.i2, self.j1:self.j2].astype(np.float)
            fld[fld==self.missing[field]] = np.nan
            setattr(self, field, fld * var.scale_factor + var.add_offset)
            return
 


    def movie(self,fieldname="salt",jd1=None,jd2=None,k=0):
        if not hasattr(self,'mp'):
            self.mp = projmaps.Projmap(self.region)
            self.xll,self.yll = self.mp(self.llon,self.llat)
            self.xl,self.yl = self.mp(
                [self.llon[0,0], self.llon[0,-1], self.llon[-1,-1],
                 self.llon[-1,0],self.llon[0,0]],
                [self.llat[0,0], self.llat[0,-1], self.llat[-1,-1],
                 self.llat[-1,0], self.llat[0,0]]
                )

        jd1 = pl.date2num(dtm(2006,9, 1))
        jd2 = pl.date2num(dtm(2006,9,30))
        for jd in np.arange(jd1,jd2,0.125):
            self.load(fieldname,jd)
            pl.clf()
            self.mp.plot(self.xl,self.yl,'0.5')
            self.mp.pcolormesh(self.xll,self.yll, miv(self.salt[k,:,:]),
                               cmap=cm.Paired)
            pl.clim(15,32)
            pl.colorbar(pad=0,aspect=40,
                        orientation="horizontal",shrink=0.5)
            print jd
            pl.title(pl.num2date(jd).strftime("%Y-%m-%d %H:%M"))
            pl.savefig('figs/movies/salt_%i.png' % int(jd*1000),dpi=100)
        
    def init_nasa(self):
        import pysea.NASA
        self.ns = pysea.NASA.nasa(res='9km')

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
        gc = pycdf.CDF(self.datadir + '/grid.cdf')
        self.landmask = gc.var('depth')[:]
        self.landmask[self.landmask>0] = 1
        self.landmask[self.landmask<0] = 0

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

    

class Sat(base.Grid):
    def __init__(self,res='250m',box=8):
        if res=='250m':
            self.satdir="/projData/casco/modis/remaps_250/box%i/" % box
        elif res=='500m':
            self.satdir="/projData/casco/modis/remaps_500/box%i/" % box
        else:
            raise "Wrong resolution"
        gridfile = glob.glob(self.satdir + "*.nav")[0]
        g = SD(gridfile, SDC.READ)
        self.llat = g.select('Latitude')[:]
        self.llon = g.select('Longitude')[:]
        self.i1=0;self.j1=0
        self.i2,self.j2 = self.llat.shape

    def read(self,filename,field='nLw_645'):
        "Read a field from a hdf file"
        h = SD(filename, SDC.READ)
        fld = h.select(field)[:]
        fld[fld<-998.] = np.nan
        self.__dict__[field] = fld
        
    def load(self,field='nLw_645',jd=0,yr=0,mn=1,dy=1):
        if jd!=0:
            yr = pl.num2date(jd).year
            mn = pl.num2date(jd).month
            dy = pl.num2date(jd).day
            yd = jd-pl.date2num(dtm(yr,1,1)) + 1
            t  = int(np.floor(np.modf(jd)[0]*8))
        elif yr!=0:
            jd = pl.date2num(yr,mn,dy)
            md  = jd - pl.date2num(dtm(1992,10,05))
        files = glob.glob(self.satdir + "A%i%03i_*.l2_remap" % (yr,yd) )
        sz = []
        for f in files:
            self.read(f,field)
            sz.append(len(
                self.__dict__[field][~np.isnan(self.__dict__[field])]))
        print sz
        self.read(files[np.nonzero(np.array(sz)==max(sz))[0][0]],field)

    def add_jds(self):
        self.jds = []
        files = glob.glob(self.satdir + "/*remap")
        for f in files:
            fp = f.split('/')[-1]
            self.jds.append(pl.date2num(
                dtm(int(fp[1:5]),1,1))+int(fp[5:8])-1)


def satmaps(res='500m',box=8, field="chlor_a"):

    cs = Sat(res=res, box=box)
    csg = GCM()
    mp = projmaps.Projmap('casco')
    x,y = mp(cs.llon,cs.llat)
    cs.xl,cs.yl = mp(
        [csg.llon[0,0], csg.llon[0,-1], csg.llon[-1,-1],
         csg.llon[-1,0],csg.llon[0,0]],
        [csg.llat[0,0], csg.llat[0,-1], csg.llat[-1,-1],
         csg.llat[-1,0], csg.llat[0,0]]
        )

    matches = glob.glob(cs.satdir + "*remap")

    for n,f in enumerate(matches):
        cs.read(f,field)
        pl.clf()
        mp.pcolormesh(x,y,miv(np.log(cs.__dict__[field])))
        mp.nice()
        pl.title(f)
        pl.clim(-3,2.5)
        mp.plot(cs.xl,cs.yl,'0.5')
        pl.savefig("/Users/bror/casfigs/cas_%03i.png" % n, dpi=100)


class GOM(base.Grid):

    def __init__(self, **kwargs):
        super(GOM, self).__init__()
        self.add_mp()
        
    def setup_grid(self):
        gc = netcdf_file(self.gridfile)
        print self.gridfile
        self.llat  = gc.variables['y'][:]
        self.llon  = gc.variables['x'][:]
        self.depth = gc.variables['depth'][:].copy()
        self.depth[self.depth<0] = np.nan
