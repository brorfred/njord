
import importlib

import pytest
import njord
import njord.nasa
from njord.utils import yrday, time
        
def test_yearday():
    tdict = njord.utils.time.expand_timeparams(jd=722815)
    assert tdict == {'yr':1980, 'mn':1, 'dy':1, 'hr':0, 'min':0, 'sec':0,
                     'jd':722815, 'yd':1.0, 'datestr':'1980-01-01'}
    ns = njord.nasa.CZCS()
    ns._timeparams(jd=722815)
    assert ns.yr == 1980
    assert ns.yd == 1

    assert (2003,1, 5) == yrday.date(yrday.jd(2003,1, 5))
    assert (2002,1, 5) == yrday.date(yrday.jd(2002,1, 5))
    assert (2004,2,28) == yrday.date(yrday.jd(2004,2,28))
    assert (2004,2,29) == yrday.date(yrday.jd(2004,2,29))
    assert (2004,3, 1) == yrday.date(yrday.jd(2004,3, 1))
    assert (2005,2,28) == yrday.date(yrday.jd(2005,2,28))
    assert (2005,3, 1) == yrday.date(yrday.jd(2005,3, 1))
    assert (2018, 12, 15, 0.5) == yrday.date(yrday.jd(2018,12,15)+0.5)

def test_datestr():
    ns = njord.nasa.CZCS()
    
    tdict = time.expand_timeparams(datestr="2018-07-01 18:00")
    assert tdict == {'yr':2018, 'mn':7, 'dy':1, 'hr':18, 'min':0, 'sec':0,
                'jd':736876.75, 'yd':182.75, 'datestr':'2018-07-01 18:00'}
  
    tdict = time.expand_timeparams(datestr="2018-07-01 18:00", hr=9)
    assert tdict == {'yr':2018, 'mn':7, 'dy':1, 'hr':9, 'min':0, 'sec':0,
                'jd':736876.375, 'yd':182.375, 'datestr':'2018-07-01 09:00'}

    tdict = time.expand_timeparams(datestr="2018-07-01 18:00", yr=2017)
    assert tdict == {'yr':2017, 'mn':7, 'dy':1, 'hr':18, 'min':0, 'sec':0,
                'jd':736511.75, 'yd':182.75, 'datestr': '2017-07-01 18:00'}

    tdict1 = time.expand_timeparams(datestr="2018-07-01 18:00", mn=8)
    tdict2 = time.expand_timeparams(datestr="2018-08-01 18:00")
    assert tdict1 == tdict2

    tdict1 = time.expand_timeparams(datestr="2018-07-01 18:00", mn=8, dy=14)
    tdict2 = time.expand_timeparams(datestr="2018-08-14 18:00")
    assert tdict1 == tdict2


def test_priority():
    ns = njord.nasa.CZCS()
    ns._timeparams(jd=722815)
    assert time.expand_timeparams(ns)["jd"] == 722815
    assert time.expand_timeparams(ns, jd=736511.75)["jd"] == 736511.75
    assert time.expand_timeparams(
        ns, jd=736511.75, datestr="2017-12-31")["jd"] == 736694.0
    assert time.expand_timeparams(
        ns, jd=736511.75, datestr="2017-12-31", hr=18)["jd"] == 736694.75
