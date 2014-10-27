#!/usr/bin/env python

import os, sys, platform
from setuptools import setup

# Version number
major = 1
minor = 3
maintenance = 9

with open("README", 'r') as file:
    readme = file.read()

# TODO: Add eventual commandline scripts here:
scripts = [
    os.path.join("scripts", "cbcflow-showcase"),
    os.path.join("scripts", "cbcflow-get-data"),
    os.path.join("scripts", "cbcflow-get-demos"),
    ]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*' % os.path.split(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)

setup(name = "cbcflow",
      version = "%d.%d.%d" % (major, minor, maintenance),
      description = "cbcflow -- Navier-Stokes solver framework from the Center of Biomedical Computing",
      long_description = readme,
      author = "Oyvind Evju, Martin Sandve Alnaes, Kent-Andre Mardal",
      author_email = "cbcflow@simula.no", 
      url = 'https://bitbucket.org/simula_cbc/cbcflow',
      classifiers = [
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Education',
          'Programming Language :: Python :: 2.7',
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
                  "cbcflow.utils",
                  "cbcflow.utils.bcs",
                  "cbcflow.utils.common",
                  "cbcflow.utils.core",
                  "cbcflow.utils.schemes",
                  ],
      package_dir = {"cbcflow": "cbcflow"},
    )

