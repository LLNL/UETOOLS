

import sys
import setuptools
from distutils.core import setup
from glob import glob


install_requires = [
     'requests', 'h5pickle'
]
if sys.version_info < (3,):
    install_requires.append('tk')
    install_requires.append('pathlib2')
elif sys.version_info > (3,0) and sys.version_info < (3,5):
    install_requires.append('tk')
    install_requires.append('pathlib')
elif sys.version_info >= (3,5):
    install_requires.append('tk')

setup(
    name='UeTools',
    version='0.1.0',
    packages=['uetools','uetools.UeGui','uetools.UeCase','uetools.UePlot'],
    maintainer='Bill Meyer',
    maintainer_email='meyer8@llnl.gov',
    description="LLNL Uedge tools",
    package_data={ "uetools":["yamls/*.yaml","yamls/*.yml"],},
    license='LLNL',
    long_description=open('README.md').read(),
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3'],
    install_requires = install_requires,
)

