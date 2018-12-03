
import numpy as np
import pylab as pl



def _semivar(distmat, fldmat1, fldmat2, lag, bw):
    """ Calculate Semivariance for a single lag """
    mask = (distmat >= (lag-bw)) & (distmat <= (lag+bw))
    return np.sum((fldmat1-fldmat2)[mask]**2) / (2.0 * np.sum(mask))


 
def semivariogram(fldvec, distmat, lag1, lag2, bw):
    '''Experimental variogram for a collection of lags
    '''
    fldmat1,fldmat2 = np.meshgrid(fldvec, fldvec)
    lagvec = np.arange(lag1, lag2+bw, bw) 
    sv = [_semivar(distmat, fldmat1, fldmat2, lag, bw) for lag in lagvec]
    return lagvec, sv
    
def C( P, h, bw ):
    '''
    Calculate the sill
    '''
    c0 = np.var( P[:,2] )
    if h == 0:
        return c0
    return c0 - SVh( P, h, bw )
