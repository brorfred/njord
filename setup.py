"""Setup file to generate a distribution of njord

usage:    python setup.py sdist
          python setup.py install
"""
from setuptools import setup, find_packages

setup(name = 'njord',
      version = '0.5.9.5',
      description = 'Package to simplify working with gridded data',
      long_description = "README.md",
      #long_description=open('docs/README.rst', 'rt').read()
      author = 'Bror Jonsson',
      author_email = 'brorfred@gmail.com',
      url = 'http://brorfred.org/python_dists/njord/',
      download_url='https://github.com/brorfred/njord/tarball/0.5',
      install_requires = ["numpy>=1.5",
                          "scipy>=0.12",
                          "matplotlib>=1.1.0",
                          "requests>=1.1.0",
                          "pyresample>=1.1.3",
                          "Pydap>=3.1",
                          "bs4>=0.0.1",
                          "netCDF4>=1.2",
                          "projmap>=0.8",
                          "click>=6.0",
                          "keyring>=13",
                          "six>1"],
      packages = ['njord','njord.utils'],
      scripts = ["njord/njordstatus",],
      package_data = {'njord': ['njord.cfg']},
     )
