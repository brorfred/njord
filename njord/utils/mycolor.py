
import pylab as pl
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize

def old_freecbar(rect,labelvec,cmap=cm.jet,label_pos='bottom',norm=None, **kwargs):
    """ Freefloating colorbar defined by rect=[lower x,lower y,len x, len y]
    Four or five labels.
    """
    numtick =  len(labelvec)
    if numtick == 6:
        tickvec = np.linspace(0,100,6)
    elif numtick == 5:
        tickvec = [0,25,50,75,100]
    elif numtick == 4:
        #tickvec = [20,40,60,80]
        tickvec = [0,33,67,100]
    elif numtick == 3:
        tickvec = [0,50,100]
    if label_pos=='left' or label_pos=='right': 
        cv = np.outer(np.ones(10),np.arange(0,1,0.01))
    else:
        cv = np.outer(np.arange(0,1,0.01),np.ones(10))
    cax = pl.axes(rect)
    if norm is None: norm = Normalize()
    pl.imshow(cv.transpose(),aspect='auto',origin='lower',cmap=cmap, norm=norm)
    if label_pos=='left' or label_pos=='right': 
        cax.yaxis.set_ticks_position(label_pos)
        cax.set_xticks([])
        cax.set_yticks(tickvec)
        cax.set_yticklabels(labelvec, **kwargs)
    else:
        cax.xaxis.set_label_position(label_pos)    
        cax.set_yticks([])
        cax.set_xticks(tickvec)
        cax.set_xticklabels(labelvec, **kwargs)
    return cax


def freecbar(rect,labelvec,cmap=cm.jet, fig=None, label_pos='bottom',
             label="", norm=None, **kwargs):
    """ Freefloating colorbar defined by rect=[lower x,lower y,len x, len y]
    Four or five labels.
    """
    fig = pl.gcf() if fig is None else fig

    numtick =  len(labelvec)
    if numtick == 6:
        tickvec = np.linspace(0,100,6)
    elif numtick == 5:
        tickvec = [0,25,50,75,100]
    elif numtick == 4:
        #tickvec = [20,40,60,80]
        tickvec = [0,33,67,100]
    elif numtick == 3:
        tickvec = [0,50,100]
    if label_pos=='left' or label_pos=='right': 
        cv = np.outer(np.ones(10),np.arange(0,1,0.01))
    else:
        cv = np.outer(np.arange(0,1,0.01),np.ones(10))

    vmin = min(labelvec)
    vmax = max(labelvec)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cax = fig.add_axes(rect)
    
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
    cb1.set_label(label)
    cb1.set_ticks(labelvec)
    cb1.ax.tick_params(**kwargs)
    pl.draw()
    return cb1

        
    if norm is None: norm = Normalize()
    pl.imshow(cv.transpose(),aspect='auto',origin='lower',cmap=cmap, norm=norm)
    if label_pos=='left' or label_pos=='right': 
        cax.yaxis.set_ticks_position(label_pos)
        cax.set_xticks([])
        cax.set_yticks(tickvec)
        cax.set_yticklabels(labelvec, **kwargs)
    else:
        cax.xaxis.set_label_position(label_pos)    
        cax.set_yticks([])
        cax.set_xticks(tickvec)
        cax.set_xticklabels(labelvec, **kwargs)
    return cax


#http://stackoverflow.com/questions/2328339/how-to-generate-n-different-colors-for-any-natural-number-n
unique_colors = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", 
                 "#FF4A46", "#008941", "#006FA6", "#A30059",
                 "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC",
                 "#B79762", "#004D43", "#8FB0FF", "#997D87",
                 "#5A0007", "#809693", "#FEFFE6", "#1B4400",
                 "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                 "#61615A", "#BA0900", "#6B7900", "#00C2A0",
                 "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                 "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                 "#300018", "#0AA6D8", "#013349", "#00846F",
                 "#372101", "#FFB500", "#C2FFED", "#A079BF",
                 "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                 "#00489C", "#6F0062", "#0CBD66", "#EEC3FF",
                 "#456D75", "#B77B68", "#7A87A1", "#788D66",
                 "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
                 "#BEC459", "#456648", "#0086ED", "#886F4C",
                 "#34362D", "#B4A8BD", "#00A6AA", "#452C2C",
                 "#636375", "#A3C8C9", "#FF913F", "#938A81",
                 "#575329", "#00FECF", "#B05B6F", "#8CD0FF",
                 "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                 "#7900D7", "#A77500", "#6367A9", "#A05837",
                 "#6B002C", "#772600", "#D790FF", "#9B9700",
                 "#549E79", "#FFF69F", "#201625", "#72418F",
                 "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                 "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
                 "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
                 "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                 "#C8D0F6", "#A3A489", "#806C66", "#222800",
                 "#BF5650", "#E83000", "#66796D", "#DA007C",
                 "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
                 "#C895C5", "#320033", "#FF6832", "#66E1D3",
                 "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]


def reverse_colourmap(cdict):
    for channel in cdict:
        arr = np.array(cdict[channel])
        arr[:,1:] = arr[::-1,1:]
        cdict[channel] = [tuple(line) for line in arr]
    return cdict

_nipy_redend_cdict = {
        'red': [(0.0, 0.0, 0.0), (0.05, 0.4667, 0.4667),
                (0.10, 0.5333, 0.5333), (0.15, 0.0, 0.0),
                (0.20, 0.0, 0.0), (0.25, 0.0, 0.0),
                (0.30, 0.0, 0.0), (0.35, 0.0, 0.0),
                (0.40, 0.0, 0.0), (0.45, 0.0, 0.0),
                (0.50, 0.0, 0.0), (0.55, 0.0, 0.0),
                (0.60, 0.0, 0.0), (0.65, 0.7333, 0.7333),
                (0.70, 0.9333, 0.9333), (0.75, 1.0, 1.0),
                (0.80, 1.0, 1.0), (0.85, 1.0, 1.0),
                (0.90, 0.8667, 0.8667), (0.95, 0.60, 0.60),
                (1.0, 0.40, 0.40)],
        'green': [(0.0, 0.0, 0.0), (0.05, 0.0, 0.0),
                (0.10, 0.0, 0.0), (0.15, 0.0, 0.0),
                (0.20, 0.0, 0.0), (0.25, 0.4667, 0.4667),
                (0.30, 0.6000, 0.6000), (0.35, 0.6667, 0.6667),
                (0.40, 0.6667, 0.6667), (0.45, 0.6000, 0.6000),
                (0.50, 0.7333, 0.7333), (0.55, 0.8667, 0.8667),
                (0.60, 1.0, 1.0), (0.65, 1.0, 1.0),
                (0.70, 0.9333, 0.9333), (0.75, 0.8000, 0.8000),
                (0.80, 0.6000, 0.6000), (0.85, 0.0, 0.0),
                (0.90, 0.0, 0.0), (0.95, 0.0, 0.0),
                (1.0, 0.0, 0.0)],
        'blue': [(0.0, 0.0, 0.0), (0.05, 0.5333, 0.5333),
                (0.10, 0.6000, 0.6000), (0.15, 0.6667, 0.6667),
                (0.20, 0.8667, 0.8667), (0.25, 0.8667, 0.8667),
                (0.30, 0.8667, 0.8667), (0.35, 0.6667, 0.6667),
                (0.40, 0.5333, 0.5333), (0.45, 0.0, 0.0),
                (0.5, 0.0, 0.0), (0.55, 0.0, 0.0),
                (0.60, 0.0, 0.0), (0.65, 0.0, 0.0),
                (0.70, 0.0, 0.0), (0.75, 0.0, 0.0),
                (0.80, 0.0, 0.0), (0.85, 0.0, 0.0),
                (0.90, 0.0, 0.0), (0.95, 0.0, 0.0),
                (1.0, 0.0, 0.0)],
        }
    
def nipy_redend():
    return mpl.colors.LinearSegmentedColormap('nipy_redend',
                                              _nipy_redend_cdict, 256)


def nipy_redend_r():
    return mpl.colors.LinearSegmentedColormap('nipy_redend_r',
                                    reverse_colourmap(_nipy_redend_cdict), 256)




def WRY():
    cdict = {'red': (  (0.00, 1.0, 1.0),
                       (0.20, 1.0, 1.0),
                       (0.60, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.20, 0.5, 0.5),
                       (0.60, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 1.0, 1.0),
                       (0.20, 0.4, 0.4),
                       (0.60, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('WhiteRedYellow',cdict,256)

def GBW():
    cdict = {'red': (  (0.00, 0.6, 0.6),
                       (0.40, 0.0, 0.0),
                       (0.80, 0.4, 0.4),
                       (1.00, 1.0, 1.0)),

             'green':( (0.00, 1.0, 1.0),
                       (0.40, 0.0, 0.0),
                       (0.80, 1.0, 1.0),
                       (1.00, 1.0, 1.0)),

             'blue': ( (0.00, 0.6, 0.6),
                       (0.40, 0.4, 0.4),
                       (0.80, 1.0, 1.0),
                       (1.00, 1.0, 1.0))}
    return mpl.colors.LinearSegmentedColormap('GreenBlueWhite',cdict,256)

def GBW_R():
    cdict = {'red': (  (0.00, 1.0, 1.0),
                       (0.20, 0.4, 0.4),
                       (0.60, 0.0, 0.0),
                       (1.00, 0.6, 0.6)),
             'green':( (0.00, 1.0, 1.0),
                       (0.20, 1.0, 1.0),
                       (0.60, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue': ( (0.00, 1.0, 1.0),
                       (0.20, 1.0, 1.0),
                       (0.60, 0.4, 0.4),
                       (1.00, 0.6, 0.6))}
    return mpl.colors.LinearSegmentedColormap('GreenBlueWhite',cdict,256)



def GBRY():
    cdict = {'red': (  (0.00, 0.6, 0.6),
                       (0.20, 0.0, 0.0),
                       (0.40, 0.4, 0.4),
                       (0.50, 1.0, 1.0),
                       (0.60, 1.0, 1.0),
                       (0.80, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.20, 0.0, 0.0),
                       (0.40, 1.0, 1.0),
                       (0.50, 1.0, 1.0),
                       (0.60, 0.5, 0.5),
                       (0.80, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 0.6, 0.6),
                       (0.20, 0.4, 0.4),
                       (0.40, 1.0, 1.0),
                       (0.50, 1.0, 1.0),
                       (0.60, 0.4, 0.4),
                       (0.80, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('GreenBlueRedYellow',cdict,256)

def GBRY9010():
    cdict = {'red': (  (0.00, 0.6, 0.6),
                       (0.40, 0.0, 0.0),
                       (0.80, 0.4, 0.4),
                       (0.91, 1.0, 1.0),
                       (0.95, 1.0, 1.0),
                       (0.97, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.40, 0.0, 0.0),
                       (0.80, 1.0, 1.0),
                       (0.91, 1.0, 1.0),
                       (0.95, 0.5, 0.5),
                       (0.97, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 0.6, 0.6),
                       (0.40, 0.4, 0.4),
                       (0.80, 1.0, 1.0),
                       (0.91, 1.0, 1.0),
                       (0.95, 0.4, 0.4),
                       (0.97, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('GreenBlue90RedYellow10',cdict,256)

def GrGr():
    cdict = {'red': (  (0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ),
             'green': ((0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ),
             'blue':  ((0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ) }
    return mpl.colors.LinearSegmentedColormap('AllGrey',cdict,8)

