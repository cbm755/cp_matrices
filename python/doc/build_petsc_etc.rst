Build instructions
==================

You need a patched PETSc version (>=3.3), petsc4py (also a patched
version) and mpi4py. In order to install that, you also need
Cython. To run our code, you need a recent NumPy version (>=1.6),
matplotlib (>=1.0) and Mayavi (?).

Everything has to be built and ran with the same Python version. If
you used EPD (free for academic usage) to get the recent
matlplotlib/NumPy versions, you will have to run all "python
setup.py..." commands from the EPD interpreter (ie, "$
/path/to/epd/install/bin/python setup.py...")

I will keep all the code under the ~/code directory, and assume that
the cp_matrices repository is also under ~/code

Cython
######

::

   curl -O http://www.cython.org/release/Cython-0.16.tar.gz
   tar xvf Cython-0.16.tar.gz
   cd Cython-0.16
   python setup.py install --user  # I use epd python (2.7.3)

PETSc
#####

::

   cd
   mkdir code/petsc
   cd code/petsc/
   hg clone http://petsc.cs.iit.edu/petsc/releases/petsc-3.3
   hg clone http://petsc.cs.iit.edu/petsc/releases/BuildSystem-3.3 petsc-3.3/config/BuildSystem
   cd petsc-3.3/
   export PETSC_ARCH=linux64
   export PETSC_DIR="$HOME/code/petsc-3.3/"
   # Patch petsc
   patch src/dm/ao/impls/mapping/aomapping.c ../cp_matrices/python/cp/petsc/patch_for_petsc.patch
   ./configure --with-shared-libraries --with-64-bit-indices --with-debugging=0 # In my laptop also: --with-c2html=0
   make PETSC_DIR=$HOME/code/petsc-3.3/ PETSC_ARCH=linux64 all
   make PETSC_DIR=$HOME/code/petsc-3.3/ PETSC_ARCH=linux64 test

PETSc4py
########

::

   cd petsc/
   hg clone https://code.google.com/p/petsc4py 
   cd petsc4py
   # Put the petsc patch in this directory, and apply it:
   patch -p1 <petsc4py_patch.patch
   python setup.py build  # again, using always the same python version
   python setup.py install --user

mpi4py
######

::

   hg clone https://code.google.com/p/mpi4py
   cd mpi4py
   python setup.py build
   python setup.py install --user


NumPy
#####

Probably the best (easiest) way to go is to get EPD. Or a newer Ubuntu distro.

If you choose to build your own numpy, do that before building petsc4py.

Thats's it!
