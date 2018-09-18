from __future__ import absolute_import, division, print_function
import importlib

import pytest
import njord

llbox = {"lat1":20, "lat2":40, "lon1":-70, "lon2":-50}


def grid_attributes(obj):
    assert obj.landmask.dtype == bool
    for attr in ["latvec","lonvec","llon","llat"]:
        assert hasattr(obj, attr)
    assert obj.shape == obj.landmask.shape
    assert obj.shape == obj.llon.shape
    assert obj.shape == obj.llat.shape
    if not hasattr(obj, "_partial_grid"):
        assert obj.landmask[0,0] == True
        assert obj.landmask[obj.jmt//2, obj.imt//2] == False
    
    assert obj.llon.shape[1] == obj.i2 - obj.i1
    assert obj.llat.shape[0] == obj.j2 - obj.j1

def config_attributes(obj):
    for attr in ["defaultjd", "minjd","maxjd","fieldlist","map_region"]:
        assert hasattr(obj, attr)

def setup_object(modulename, classname, **kwargs):
    clss = getattr(importlib.import_module("njord.%s" % modulename), classname)
    obj = clss(**kwargs)
    obj.datadir = "/tmp/njord/tests/"
    config_attributes(obj)
    grid_attributes(obj)
    for fieldname in obj.fieldlist:
        fld = obj.get_field(fieldname)
        assert fld.shape == obj.shape
    return obj
    
    
def test_nasa_9km():
    for instrument in ["MODIS", "VIIRS", "SeaWiFS", "OCTS", "CZCS"]:
        ns = setup_object("nasa", instrument)
        ns = setup_object("nasa", instrument, **llbox)
        assert ns.shape == (242, 242)

def test_nasa_4km():
    for instrument in ["MODIS", "VIIRS"]:
        ns =  setup_object("nasa", instrument, res="4km")
        ns =  setup_object("nasa", instrument, res="4km", **llbox)
        assert ns.shape == (482, 482)

def test_seawinds():
        ns = setup_object("winds", "Seawinds")
        ns = setup_object("winds", "Seawinds", **llbox)
        assert ns.shape == (83, 81)

def test_ccmp():
        ns = setup_object("winds", "CCMP")
        ns = setup_object("winds", "CCMP", **llbox)
        assert ns.shape == (82, 82)
        
def test_ncep():
        ns = setup_object("winds", "NCEP")
        ns = setup_object("winds", "NCEP", **llbox)
        assert ns.shape == (11, 9)

def test_L4_K10():
        ns = setup_object("ghrsst", "L4_K10")
        ns = setup_object("ghrsst", "L4_K10", **llbox)
        assert ns.shape == (200, 201)

def test_L4CMC():
        ns = setup_object("ghrsst", "L4CMC")
        ns = setup_object("ghrsst", "L4CMC", **llbox)
        assert ns.shape == (103, 101)


def test_cci():
        ns = setup_object("cci", "OceanColor")
        ns = setup_object("cci", "OceanColor", **llbox)
        assert ns.shape == (482, 482)

        for timetype in ns.timetypelist:
            ns.load(timetype=timetype)
            ns.chl.shape == (482, 482)