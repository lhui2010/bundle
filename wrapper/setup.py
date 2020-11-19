#!/usr/bin/env python

from setuptools import setup, find_packages
import os.path as op
#import versioneer


name = "iga"

setup_dir = op.abspath(op.dirname(__file__))
requirements = [x.strip() for x in open(op.join(setup_dir, "requirements.txt"))]
packages = [name] + [
    ".".join((name, x)) for x in find_packages("jcvi", exclude=["test*.py"])
]

setup(
    name = name,
    version = "0.1",
    packages = packages,
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    # metadata for upload to PyPI
    author = "Me",
    author_email = "me@example.com",
    description = "iga",
    license = "BSD",
    install_requires=requirements,
)
