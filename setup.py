#import os
import glob
#import multiprocessing # Why is this line necessary?
from setuptools import setup#, find_packages

#def read(fname):
#    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='sptspd',
    version='0.1',
    description='My current mystery project',
    #long_description=read('README.md'),
    #packages = find_packages(),
    author='Michael Gruenstaeudl',
    author_email='m.gruenstaeudl@fu-berlin.de',
    #url='https://github.com/michaelgruenstaeudl/SPTSPD',
    scripts=glob.glob('scripts/*'),
    #test_suite='nose.collector',
    #tests_require=['nose >= 1.3', 'mock'],
    license='GPLv3',
    packages=['sptspd'],
    #classifiers=[
    #    "License :: OSI Approved :: GNU General Public License (GPLv3)",
    #    "Programming Language :: Python",
    #    "Development Status :: 1 - Planning",
    #    "Intended Audience :: Science/Research",
    #    "Topic :: Scientific/Engineering :: Bio-Informatics",
    #    ],
    zip_safe=False
)
