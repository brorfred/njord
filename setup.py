"""Setup file to generate a distribution of njord

usage:    python setup.py sdist
          python setup.py install
"""


from distutils.core import setup

setup(name = 'njord',
      version = '0.5.9.1',
      description = 'Package to simplify working with gridded data',
      long_description = "README.md",
      #long_description=open('docs/README.rst', 'rt').read()

      author = 'Bror Jonsson',
      author_email = 'brorfred@gmail.com',
      url = 'http://brorfred.org/python_dists/njord/',
      requires = ["numpy(>=1.5)",
                  "matplotlib(>=1.1.0)",
                  "mpl_toolkits(>=1.0)",
                  "requests(>=1.1.0)",
                  "projmap(>=0.5)",
                  "requests(>=1.1)"],
      packages = ['njord'],
      scripts = ["njord/njordstatus",],
      package_data = {'njord': ['njord.cfg']},
     )
