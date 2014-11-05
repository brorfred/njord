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
                          "matplotlib>=1.1.0",
                          "requests>=1.1.0",
                          "projmap>=0.5",],
      packages = ['njord'],
      scripts = ["njord/njordstatus",],
      package_data = {'njord': ['njord.cfg']},
     )
