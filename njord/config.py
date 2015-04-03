import os
import json
import ConfigParser

def load(filepref, projname, kwargs={}):
    """Read and parse the config file"""
    preset_dict = {}
    basedir =  os.path.dirname(os.path.abspath(__file__))

    cfg = ConfigParser.ConfigParser()
    files = ["%s/%s.cfg" %  (os.curdir, filepref),
             "%s/.%s.cfg" % (os.path.expanduser("~"), filepref),
             "%s/%s.cfg" %  (basedir, filepref)]
    
    for fnm in files:
        cfg.read(fnm)
        if projname in cfg.sections():
            preset_dict['config_file'] = fnm
            break
    else:
        raise NameError('Project not included in config files')
        
    def splitkey(key, val):
        if key in kwargs.keys():
            preset_dict[key] = kwargs[key]
            del kwargs[key]
        else:
            preset_dict[key] = val

    for key,val in cfg.items(projname):
            if type(val) is str and "~/" in val:
                val = os.path.expanduser(val)
            try:
                splitkey(key, json.loads(val))
            except ValueError:
                splitkey(key, val)
                
    preset_dict.update(kwargs)
    return preset_dict
