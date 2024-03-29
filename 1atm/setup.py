from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
      name = "HMX",
      ext_modules = [
      Extension("melt_layer",
               sources = ["melt_layer.pyx","Phase.cpp","Species.cpp", "Reactions.cpp"],
               language = "c++"),
      ],
      cmdclass = {'build_ext':build_ext},    
)

