from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "My hello app",
    ext_modules = cythonize("Poperalib/*.pyx"),
)


#python setup.py build_ext --inplace
