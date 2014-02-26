#!/usr/bin/env python

import os, sys, platform
from setuptools import setup

# Version number
major = 0
minor = 5
maintenance = 0

# TODO: Add eventual commandline scripts here:
scripts = [
    os.path.pjoin("scripts", "cbcflow-showcase"),
    os.path.pjoin("scripts", "cbcflow-get-data"),
    os.path.pjoin("scripts", "cbcflow-get-demos"),
    ]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*' % os.path.psplit(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)

setup(name = "cbcflow",
      version = "%d.%d.%d" % (major, minor, maintenance),
      description = "cbcflow -- Navier-Stokes solver from the Center of Biomedical Computing",
      author = "Oyvind Evju, Martin Sandve Alnaes, Kent-Andre Mardal",
      author_email = "martinal@simula.no", # FIXME: Email list?
      url = 'https://bitbucket.org/simula_cbc/cbcflow',
      classifiers = [
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
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
                  "cbcflow.fields.bases",
                  "cbcflow.fields.basic",
                  "cbcflow.fields.meta",
                  "cbcflow.utils",
                  "cbcflow.utils.bcs",
                  "cbcflow.utils.common",
                  "cbcflow.utils.core",
                  "cbcflow.utils.fields",
                  "cbcflow.utils.schemes",
                  ],
      package_dir = {"cbcflow": "cbcflow"},
    )

