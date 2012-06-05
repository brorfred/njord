
import numpy as np
import pylab as pl

import projmaps


def llmat():

    """
    http://spg.ucsd.edu/Satellite_Projects/CAL/
    http://cwcaribbean.aoml.noaa.gov/hdf.html
    http://elonen.iki.fi/code/misc-notes/affine-fit/
    """
    lon = np.linspace(-135,-100,3840)
    lat = np.linspace(16,45,3405)[::-1]
    llon,llat = np.meshgrid(lon, lat)

    ivec = np.arange(3405)
    jvec = np.arange(3840)
    imat,jmat = np.meshgrid(ivec,jvec)
    
    et_affine=[-1.00034410428395, 0, 0, 1.00034410428395, 1703, 1920.5]
    a = et_affine[0] 
    b = et_affine[1]
    c = et_affine[2]
    d = et_affine[3]
    e = et_affine[4]
    f = et_affine[5]

    imatp = imat*a + jmat*b + e
    jmatp = imat*c + jmat*d + f

    


    return llon,llat
