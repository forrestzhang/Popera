# -*- coding: utf-8 -*-

import os
from setuptools import setup
from setuptools import find_packages

def readme(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='Popera',
      version='1.0.3',
      packages=find_packages(),
      author='Tao Zhang',
      author_email='zhangtao@yzu.edu.cn',
      scripts=['bin/popera', 'bin/popera_dhs_count'],
      py_modules=['Popera', 'Popera_dhs_count'],
      include_package_data=True, 
      url='https://github.com/forrestzhang/Popera',
      license='LICENSE',
      description='Popera: DNase I hypersensitive sites identification',
      long_description=readme('README.md'),
      classifiers=['Intended Audience :: Science/Research'],
      install_requires=["numpy", "scipy", "pysam", "pyBigWig"],
      zip_safe=False,
     )
