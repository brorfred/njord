"""Setup file to generate a distribution of projmap

use as python setup.py sdist"""


from distutils.core import setup

setup(name='njord',
      version='0.5',
      description='Package to simplify working with gridded data',
      author='Bror Jonsson',
      author_email='brorfred@gmail.com',
      url='http://brorfred.org/python_dists/njord/',
      requires=["numpy(>=1.5)", "matplotlib(>=1.1.0)", "mpl_toolkits(>=1.0)",
                "projmap(>=0.5)"],
      packages=['njord'],
      package_data={'njord': ['njord.cfg']},
     )
