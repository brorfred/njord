njord
=====
A package to generalize and simplify working with different types of 
geophysical/oceanographical/meterological data. The package allows for 
analysis and presentation of data from different sources in a unified
way and provides a number of functions for this. The package requires
some work to when first setup, both by editing the config file and 
most likely adding modules/classes. 

The structure is as follows: 

A class represents a specific dataset.
A module contains all classes from the same source. They are normally
setup in the same format. 


Installation
------------
 - (Edit config file and add modules/classes.)
 - sudo python setup.py install      # global installation
 - python setup.py install --user    # installation for current user


Njord can reproject beween different grids by using the *pyresample* package from pytroll.org. This package is not installed by default due to some problems with some of it's prerequisites on Macos. You can install it with the command:

    pip install pyresample


If you want to use njord with the repqojection capality in MacOS, you might have problems with Apple's installation of gcc not finding -fopenmp when trying to install pyresample. If so, try to use LLVM installed via homebrew instead:

    brew install llvm
    export CC=/usr/local/opt/llvm/bin/clang
    pip install pyresample

Define a project in one of the following files:

./njord.cfg
~/.njord.cfg
/path/to/package/njord.cfg

Example:

```ini
[Default]
basedir:     /projData

[rutgers.Coral]
datadir:     %(basedir)s/rutgers/CORAL/
gridfile:    %(datadir)s/coral_grd.nc
map_region:  indthr
imt:         1281
jmt:         641
```

Usage
-----
```python
>>> from njord import rutgers 
>>>
>>> #Create a grid instance
>>> mp = rutgers.Coral('handle')
>>>
>>> #Load u-velociites
>>> mp.load('u')
```



https://stackoverflow.com/questions/3764291/checking-network-connection