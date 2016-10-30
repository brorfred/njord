from numpy import *
import numpy as np


pi180 =  pi/180
earth_radius = 6371.0*1000
 

def lldist(lon,lat):
    
    lat = lat * pi180
    dlon = (lon[1:] - lon[:-1])*pi180 
    dlat = (lat[1:] - lat[:-1])

    a = (sin(dlat/2))**2 + cos(lat[:-1]) * cos(lat[1:]) * (sin(dlon/2))**2
    angles = lon * 0
    angles[1:] = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles

def ll2dist(lon,lat):

    lat = lat * pi180
    dlon = (lon[1,:] - lon[0,:])*pi180 
    dlat = (lat[1,:] - lat[0,:])

    a = (sin(dlat/2))**2 + cos(lat[0,:]) * cos(lat[1,:]) * (sin(dlon/2))**2
    angles = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles

def ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2):
    """Find distance between each respective elements in two sets of vectors"""
    latvec1 = np.deg2rad(latvec1)
    latvec2 = np.deg2rad(latvec2)
    lonvec1 = np.deg2rad(lonvec1)
    lonvec2 = np.deg2rad(lonvec2)

    dlon = (lonvec2 - lonvec1)
    dlat = (latvec2 - latvec1)

    a = (sin(dlat/2))**2 + cos(latvec1) * cos(latvec2) * (sin(dlon/2))**2
    angles = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles


def spherical_dist(pos1, pos2, r=3958.75):
    pos1 = np.deg2rad(pos1)
    pos2 = np.deg2rad(pos2)
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return earth_radius * np.arccos(cos_lat_d - cos_lat1 *
                                    cos_lat2 * (1 - cos_lon_d))



def distance_matrix(lonvec, latvec):
    """Calculate a distance matrix for all combinations of the lat-lon positions"""
    latmat1,latmat2 = np.meshgrid(latvec, latvec)
    lonmat1,lonmat2 = np.meshgrid(lonvec, lonvec)
    pos1 = np.vstack((latmat1.flatten(), lonmat1.flatten())).T
    pos2 = np.vstack((latmat2.flatten(), lonmat2.flatten())).T
    
    dist =  spherical_dist(pos1,pos2)
    return dist.reshape(int(np.sqrt(len(dist))),-1)
