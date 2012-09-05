from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension('cp.surfaces.triangulation_fast', 
                         ["cp/surfaces/triangulation_fast.pyx"])]

setup(
    name = "cp",
    version = "0.2dev",
    packages = ["cp",
                "cp.tools",
                "cp.petsc",
                "cp.surfaces",
                "cp.tests",
                ],
    package_data = {'cp.tests': ['data/*']},
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    requires=['numpy (>=1.6)', 'scipy', 'matplotlib (>=1.0)', 'cython (>=0.16)']
)
