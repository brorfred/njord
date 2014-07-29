"""Convert a global grid where the center isnt at GMT to one."""
import numpy as np


class Shift:

    def __init__(self, lonvec, idim=1):
        """Find the new centerpoint and adjust the lon vector""" 
        self.imt = len(lonvec)
        lonvec[lonvec>180] = lonvec[lonvec>180]-360
        lonvec[lonvec<-180] = lonvec[lonvec<-180]+360
        try:
            gr = np.nonzero( (lonvec[1:]-lonvec[:-1]) < 0)[0].item() + 1
        except ValueError:
            raise GridError
        self.lonvec = np.roll(lonvec,self.imt - gr)
        self.idim = idim
        self.gmtpos = gr
    
    def field(self, fld):
        """Convert a field to GMT grid. """
        axis = np.nonzero(np.array(fld.shape) == len(self.lonvec))[0][0]
        if self.idim == 1:
            fld = np.squeeze(fld)
            if fld.shape[axis] == self.imt:
                return np.roll(fld,self.imt-self.gmtpos,axis=axis)
            else:
                raise ValueError, "Lonvec and iaxis of field not of same length"

class GridError(Exception):
    """Exception raised if Shift can't find a minimum"""
    def __init__(self):
        self.value = "Couldn't shift the lon- vecctor, center already at GMT?"
    def __str__(self):        
        return repr(self.value)
