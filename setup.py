#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""
[License: GNU General Public License v3 (GPLv3)]
"""

import blc
from setuptools import setup


def get_requirements():
    with open('requirements.txt', 'r') as fh:
        content = fh.read().strip().split()
    return content


setup(name="bam-lorenz-coverage",
      scripts=['bin/bam-lorenz-coverage'],
      packages=["blc"],
      test_suite="tests",
      tests_require=['nose2', 'pytest', 'pytest-cov', 'flake8'],
      setup_requires=[get_requirements()],
      install_requires=[get_requirements()],
      version=blc.__version__,
      description="BAM coverage stats plots",
      author=blc.__author__,
      url=blc.__homepage__,
      keywords=["NGS", "BAM", "Lorenz", "coverage"],
      classifiers=['Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ])
