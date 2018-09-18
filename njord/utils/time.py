
from datetime import datetime as dtm

from . import mpl_dates as pl 

def expand_timeparams(obj=None, datestr=None, **kwargs):
    """Expand all time parameters from given values"""

    if (obj is not None) and not hasattr(obj,"projname"):
        raise TypeError("Obj must be a njord class instance")
    if (datestr is not None) and (type(datestr)!=str):
        raise TypeError("datestr must be a date string")

    def expand_jd():
        for key1,key2 in zip(tlist,
                             ["year","month","day","hour","minute", "second"]):
            tdict[key1] = getattr(pl.num2date(tdict["jd"]), key2)
        
    tlist = ["yr", "mn", "dy", "hr", "min", "sec", "jd", "yd"]
    tdict = {key:0 for key in tlist}
    if obj is not None:
        for key in tlist + ["defaultjd"]:
            tdict[key] = getattr(obj, key, 0)
    if tdict["jd"] is 0:
        tdict["jd"] = tdict.get("defaultjd")
    if kwargs.get("jd") is not None:
        tdict["jd"] = kwargs["jd"]
    if datestr is not None:
        tdict["jd"] = pl.datestr2num(datestr)
    if tdict.get("jd") is not 0:
        expand_jd()
    for key in kwargs:
        tdict[key] = kwargs.get(key)
    tdict["jd"] = pl.date2num(
        dtm(*[tdict[key] for key in ["yr","mn","dy","hr","min","sec"]]))
    if ('yd' in kwargs):
        if tdict["yd"] < 1:
            tdict["yr"] -= 1
            ydmax = (pl.date2num(dtm(tdict["yr"], 12, 31)) -
                     pl.date2num(dtm(tdict["yr"],  1,  1))) + 1    
            tdict["yd"]  = ydmax + tdict["yd"]
        tdict["jd"] = tdict["yd"] + pl.date2num(dtm(tdict["yr"],1,1)) - 1
    else:
        tdict["yd"] = tdict["jd"] - pl.date2num(dtm(tdict["yr"],1,1)) + 1
    expand_jd()

    if tdict["jd"]==int(tdict["jd"]):
        tdict["jd"] = int(tdict["jd"])
        tdict["datestr"] = pl.num2date(tdict["jd"]).strftime("%Y-%m-%d")
    else:
        tdict["datestr"] = pl.num2date(tdict["jd"]).strftime("%Y-%m-%d %H:%M")
    return tdict

    """
    if hasattr(params,'hourlist'):
        dd = params.jd-int(params.jd)
        ddlist = np.array(params.hourlist).astype(float)/24
        ddpos = np.argmin(np.abs(ddlist-dd))
        params.jd = int(params.jd) + ddlist[ddpos]

        
    if params.jd < params.fulltvec.min():
        raise ValueError("Date before first available model date")
    if params.jd > params.fulltvec.max():
        raise ValueError("Date after last available model date")
    """
