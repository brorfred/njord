
import numpy as np
import pylab as pl

bore_ydmn    = np.cumsum(np.array([0,31,28,31,30,31,30,31,31,30,31,30]))
astr_ydmn =  np.cumsum(np.array([0,31,31,30,31,30,31,28,31,30,31,30]))-182.5

ydmnstr = ['Jan','Feb','Mar','Apr','May','Jun',
           'Jul','Aug','Sep','Oct','Nov','Dec']

mnstr = ydmnstr

def calc_parts(jd):
    l = jd.astype(np.int)+68569 + 1721425
    n = 4*l/146097
    l = l-(146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l-1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l
    return i,j,k,l

def date(jd):
    jd = np.array(jd)
    part = jd - jd.astype(np.int)
    l = jd.astype(np.int)+68569 + 1721425
    n = 4*l/146097
    l = l-(146097*n+3)/4
    yr = 4000*(l+1)/1461001
    l = l-1461*yr/4+31
    mn = 80*l/2447
    dy = l-2447*mn/80
    l = mn/11
    mn = mn+2-12*l
    yr = 100*(n-49)+yr+l
    return yr,mn,dy,part

def jd(yr, mn, dy):
    return (dy - 32075 + 1461 * (yr + 4800 + (mn - 14) / 12) / 4 +
            367 * (mn - 2 - (mn - 14) / 12 * 12) / 12 -
            3 * ((yr + 4900 + (mn - 14) / 12) / 100) / 4) - 1721427.0

def yd(jd):
    """Calculate year day from julian date. """
    l = jd.astype(np.int)+68569 + 1721425
    n = 4*l/146097
    l = l-(146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l-1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l
    j = j * 0 + 1
    k = k * 0 + 1
    jdyr = (k- 32075 + 1461 * (i+4800 + (j-14) / 12) / 4 +
            367 * (j-2-(j-14)/12*12)/12 -
              3 * ((i+4900+(j-14)/12) / 100) / 4) - 1721428
    return jd - jdyr

def asyd(jd):
    ast_yd = yd(jd)
    ast_yd[ast_yd>182.5] = ast_yd[ast_yd>182.5] - 365
    return ast_yd


def mndy(yr,yd):
    """Calculate month and day from year day. """
    jdtot = jd(yr,1,1) + yd 
    yr,mn,dy,_ = date(jdtot)
    return mn,dy


def years(jd):
    """Calculate years from a vector with julian dates. """
    i,j,k,l = calc_parts(jd)
    return i

def months(jd):
    """Calculate months from a vector with julian dates. """
    i,j,k,l = calc_parts(jd)
    return j

def xaxis():
    pl.xlim(0,365)
    pl.xticks(bore_ydmn, ydmnstr)
    pl.xlabel('Climatological Year')


def astr_xaxis(skip=1):
    pl.xlim(-182.5,182.5)
    pl.xticks(astr_ydmn[::skip], ydmnstr[6::skip] + ydmnstr[:6:skip])
    pl.xlabel('Climatological Year')

def djdticks(maxdjd):
    """Generate ticks that scale with days,weeks,years """
    maxdjd = float(maxdjd)
    if maxdjd <= 60:
        maxval = np.ceil((maxdjd/7)/2)*2
        return ['%gw' % f for f in np.linspace(0,maxval,5)], maxval*7
    if maxdjd <= 365:
        maxval = np.ceil((maxdjd/30.5)/2)*2
        return ['%gm' % f for f in np.linspace(0,maxval,5)], maxval*30.5
    maxval = np.ceil((maxdjd/365)/2)*2
    return ['%gyr' % f for f in np.linspace(0,maxval,5)], maxval*365
