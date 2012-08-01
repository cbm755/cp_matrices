from distutils.core import setup

setup(
    name="cp",
    version="0.2dev",
    packages=["cp",
              "cp.tools",
              "cp.petsc",
              "cp.surfaces",
              "cp.tests",],
    package_data={'cp.tests': ['data/*']},
)
