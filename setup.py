import os
import glob
import unittest
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='*_test.py')
    return test_suite

setup(
    name='annonex2embl',
    version='0.2',
    description='Converts annotated DNA sequence alignments in NEXUS format to EMBL submission files',
    long_description=read('README.rst'),
    classifiers=[
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: GNU General Public License (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    keywords='EMBL DNA sequence submission',
    url='https://github.com/michaelgruenstaeudl/SPTSPD',
    author='Michael Gruenstaeudl',
    author_email='m.gruenstaeudl@fu-berlin.de',
    license='GPLv3',
    packages=['annonex2embl'], # So that the subfolder 'annonex2embl' is read immediately.
    #packages = find_packages(),
    install_requires=['biopython'],
    scripts=glob.glob('scripts/*'),
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False
)
