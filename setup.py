#!/usr/bin/env python
from setuptools import setup, find_packages
import mistgen

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()
    
setup(
    name="mistgen",
    version=mistgen.__version__,
    author="Zhaoliang Zheng",
    author_email="zhz03@g.ucla.edu",
    description="minimum snap trajectory generator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GNU GENERAL PUBLIC LICENSE V3",
    url="https://github.com/zhz03/MiSTGen",
    packages=find_packages(),
    install_requires=[
        "matplotlib==3.4.2",
        "numpy",
        "cvxopt"
        ],
    classifiers=[
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
    ],
)