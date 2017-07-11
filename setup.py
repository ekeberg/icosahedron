from distutils.core import setup, Extension
import numpy.distutils.misc_util

ext = Extension("icosahedron", sources=["icosahedronmodule.c"], include_dirs=[numpy.get_include()])
setup(name="icosahedron", ext_modules=[ext], include_dirs=numpy.get_include())
