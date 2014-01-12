"""Convert a global grid where the center isnt at GMT to one."""
import numpy as np


class Shift:

    def __init__(self, lonvec, idim=1):
        """Find the new centerpoint and adjust the lon vector""" 
        lonvec[lonvec>180] = lonvec[lonvec>180]-360
        lonvec[lonvec<-180] = lonvec[lonvec<-180]+360
        try:
            gr = np.nonzero( (lonvec[1:]-lonvec[:-1]) < 0)[0].item() + 1
        except ValueError:
            raise GridError
        tmp = lonvec[:gr].copy()
        self.lonvec = np.hstack( (lonvec[gr:],tmp) )
        self.idim = idim
        self.gmtpos = gr
    
    def field(self, fld):
        """Convert a field to GMT grid. """
        if self.idim == 1:
            fld = np.squeeze(fld)
            tmp = np.squeeze(fld[...,:self.gmtpos].copy())
            return np.concatenate( (np.squeeze(fld[..., self.gmtpos:]), tmp), 
                                   axis=fld.ndim-1 )

class GridError(Exception):
    """Exception raised if Shift can't find a minimum"""
    def __init__(self):
        self.value = "Couldn't shift the lon- vecctor, center already at GMT?"
    def __str__(self):        
        return repr(self.value)
