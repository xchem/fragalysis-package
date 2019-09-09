# Fragalysis

[![Build Status](https://travis-ci.org/xchem/fragalysis.svg?branch=master)](https://travis-ci.org/xchem/fragalysis)
[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.0.38-blue.svg?style=flat)](https://github.com/xchem/fragalysis)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis/blob/master/LICENSE.txt)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis/context:python)

![PyPI](https://img.shields.io/pypi/v/fragalysis)

Basic RDKit based Python tools for analysis of protein-ligand interactions.

Currently contains: -

1.  Clustering - based on WONKA method - but separated from that code-base.
    Cluster waters, residues, ligands and pharmacophores. (Under development)
2.  Astex Fragment Network - implementation on the basis of their recent paper
3.  Conformer generation code - based on known X-ray structures

## Publishing (to PyPi)
Armed with PyPi account credentials, ideally from within a Python 3
virtual environment, run the following from the project root: -

    $ export TWINE_USERNAME=PyPiUsername
    $ export TWINE_PASSWORD=PyPiPassword
    $ pip install -r package-requirements.txt
    $ rm dist/*
    $ python setup.py bdist_wheel
    $ twine upload dist/*
