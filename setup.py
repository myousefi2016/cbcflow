#!/usr/bin/env python
#from distutils.core import setup
#from distutils import sysconfig
from setuptools import setup, find_packages
from os.path import join as pjoin, split as psplit
from glob import glob
import sys
import platform
from setuptools.command.install import install
import shutil
import subprocess

# Version number
major = 0
minor = 5
maintenance = 0

"""
class MyInstall(install):

    def run(self):
        # Install
        install.run(self)

        # Clean up after install
        subprocess.call("python setup.py clean --all", shell=True)
        shutil.rmtree("cbcflow.egg-info")
"""
        

# TODO: Add eventual commandline scripts here:
scripts = [
    pjoin("scripts", "cbcflow-showcase"),
    pjoin("scripts", "cbcflow-get-data"),
    pjoin("scripts", "cbcflow-get-demos"),
    ]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*' % psplit(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)

setup(name = "cbcflow",
      version = "%d.%d.%d" % (major, minor, maintenance),
      description = "cbcflow -- Navier-Stokes solver from the Center of Biomedical Computing",
      author = "Oyvind Evju, Martin Sandve Alnaes, Kent-Andre Mardal",
      author_email = "martinal@simula.no", # TODO: Email list?
      url = 'https://bitbucket.org/simula_cbc/cbcflow',
      classifiers = [
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2.6',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
      scripts = scripts,
      packages = ["cbcflow",
                  "cbcflow.core",
                  "cbcflow.schemes",
                  "cbcflow.schemes.official",
                  "cbcflow.schemes.experimental",
                  "cbcflow.bcs",
                  "cbcflow.fields",
                  "cbcflow.fields.bases",
                  "cbcflow.fields.basic",
                  "cbcflow.fields.meta",
                  "cbcflow.utils",
                  "cbcflow.utils.bcs",
                  "cbcflow.utils.common",
                  "cbcflow.utils.core",
                  "cbcflow.utils.fields",
                  "cbcflow.utils.schemes",
                  #"cbcflow.utils.fenicstools",
                  # Note: These are not python packages:
                  #"cbcflow.utils.fenicstools.fem",
                  #"cbcflow.utils.fenicstools.Probe",
                  ],
      package_dir = {"cbcflow": "cbcflow"},
      
      # Require and fetch fenicstools
      install_requires = ['fenicstools<=1.3.0'],
      dependency_links=['git+https://github.com/oyvinev/fenicstools.git@1.3.0-mod#egg=fenicstools'],
      
      # Use custom install class for proper cleanup
      # FIXME: This has to be handled differently, possibly makefile?
      #cmdclass={'install': MyInstall},

    )

