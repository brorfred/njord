
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

import requests
import click


def isncfile(filename):
    if not os.path.isfile(filename):
        return False
    try:
        netCDF4.Dataset(filename)
        return True
    except:
        os.unlink(filename)
        return False

def retrive_ftp(url, local_dir="./", username="anonymous", password="njord@bror.us",
                verbose=False):
    """Retrive file from remote server via http or ftp

    https://stackoverflow.com/questions/51684008
    """
    def vprint(text):
        if verbose:
            print(text)
    
    spliturl = urlsplit(url)
    vprint(spliturl)
    host = spliturl.netloc
    path = spliturl.path
    ftpdir,lfn = os.path.split(path)
    local_filename = os.path.join(local_dir, lfn)

    #if spliturl.scheme == "ftp":
    try:
        ftp = ftplib.FTP(host) 
        vprint(ftp.login(username, password))
        vprint("Change dir to '%s'" % ftpdir)
        vprint(ftp.cwd(ftpdir))
    except ftplib.error_perm as err:
        print (spliturl.netloc)
        print (os.path.split(spliturl.path)[0])
        raise IOError(err)

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

def retrive_http(url, local_filename=None, username=None, params=None, verbose=False):
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

