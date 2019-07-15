import setuptools
import sys
from setuptools import setup, find_packages, Extension

# with open("README.md", "r") as fh:
#     long_description = fh.read()
    
required = ["glypy","pandas","numpy","seaborn","networkx","ndex","xlrd"]

setuptools.setup(
     name='glycompare',  
     version='0.3.2',
     install_requires=required,
     author="B.B., B.P.K",
     author_email="bobao@eng.ucsd.edu",
     description="GlyCompare",
     long_description='Please check https://github.com/LewisLabUCSD/GlyCompare',
   long_description_content_type="text/markdown",
     url="https://github.com/LewisLabUCSD/GlyCompare",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )
