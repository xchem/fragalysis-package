# Fragalysis Package

[![build latest](https://github.com/xchem/fragalysis-package/actions/workflows/build-latest.yml/badge.svg)](https://github.com/xchem/fragalysis-package/actions/workflows/build-latest.yml)
[![pip release](https://github.com/xchem/fragalysis-package/actions/workflows/pip-release.yml/badge.svg)](https://github.com/xchem/fragalysis-package/actions/workflows/pip-release.yml)

[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-package/blob/master/LICENSE.txt)
[![CodeFactor](https://www.codefactor.io/repository/github/xchem/fragalysis-package/badge)](https://www.codefactor.io/repository/github/xchem/fragalysis-package)

![PyPI](https://img.shields.io/pypi/v/fragalysis)

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Basic RDKit based Python tools for analysis of protein-ligand interactions.

Currently contains: -

1. Clustering - based on WONKA method - but separated from that code-base.
    Cluster waters, residues, ligands and pharmacophores. (Under development)
2. Astex Fragment Network - implementation on the basis of their recent paper
3. Conformer generation code - based on known X-ray structures
4. Support for the neo4j 4.4.2 graph database

## Pre-commit
The project uses [pre-commit] to enforce linting of files prior to committing
them to the upstream repository.

To get started review the pre-commit utility and then set-up your local clone
by following the **Installation** and **Quick Start** sections of the
pre-commit documentation.

Ideally from a Python environment...

    python -m venv venv
    source venv/bin/activate

    pip install --upgrade pip
    pip install -r build-requirements.txt
    pre-commit install -t commit-msg -t pre-commit

Now the project's rules will run on every commit and you can check the
state of the repository as it stands with...

    pre-commit run --all-files

## Publishing (to PyPI)
If the repository has been provided with a `PYPI_APIKEY` **Secret**
the GitHub `pip-release` workflow action will automatically publish the package to
PyPI when you create a new **Release**. The package version will be set using
the release `tag_name`.

---

[pre-commit]: https://pre-commit.com
