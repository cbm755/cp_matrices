from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('cp.surfaces.triangulation_fast', 
                         ["cp/surfaces/triangulation_fast.pyx"])]

setup(
    name="cp",
    version="0.2dev",
    packages=["cp",
              "cp.tools",
              "cp.petsc",
              "cp.surfaces",
              "cp.tests",],
    package_data={'cp.tests': ['data/*']},
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
