Build instructions
==================

You need a patched PETSc version (>=3.3), petsc4py (also a patched
version) and mpi4py. You also need Cython, a recent NumPy version
(>=1.6), matplotlib (>=1.0), Mayavi and scipy.

Everything has to be built and ran with the same Python version. If
you are using EPD (free for academic usage) in your workstation or
cyclops to get the recent matlplotlib/NumPy versions, make sure you're
always using EPD python.

Be careful with you PYTHONPATH. Most of the trouble we had getting the
code to run in Yujia's workstation was due to a non empty PYTHONPATH:
we were picking up old packages from SAGE.

I will keep all the code under the ~/code directory, and assume that
the cp_matrices repository is also under ~/code

Cython
######

Cython 0.16 is in EPD 7.3, and that should be enough.

Else, and if pip is not available in your machine, you can always do

::

   curl -O http://www.cython.org/release/Cython-0.17.tar.gz
   tar xvf Cython-0.17.tar.gz
   cd Cython-0.17
   python setup.py install --user

Documentation available in http://docs.cython.org/src/quickstart/install.html

PETSc
#####

::

   cd
   mkdir code/petsc
   cd code/petsc/
   hg clone http://petsc.cs.iit.edu/petsc/releases/petsc-3.3  # Takes a while
   hg clone http://petsc.cs.iit.edu/petsc/releases/BuildSystem-3.3 petsc-3.3/config/BuildSystem
   cd petsc-3.3/
   export PETSC_ARCH=linux64  # or something else
   export PETSC_DIR="$HOME/code/petsc-3.3/"  # Adjust accordingly
   patch src/dm/ao/impls/mapping/aomapping.c ../cp_matrices/python/patches/petsc.patc
   patch src/dm/ao/impls/mapping/aomapping.c ../cp_matrices/python/patches/petsc_2.patch

For the configure step, you can disable generating documentation using
``--with-c2html=0``. ``--with-64-bit-indices`` is needed in 64 bit
OSes as far as I know (ie, I need it in my laptop, workstation and
cyclops)::

   ./configure --with-shared-libraries --with-64-bit-indices --with-debugging=1
   make PETSC_DIR=$HOME/code/petsc-3.3/ PETSC_ARCH=linux64 all
   make PETSC_DIR=$HOME/code/petsc-3.3/ PETSC_ARCH=linux64 test

You can add this to your `.bashrc`::

   export PETSC_DIR=$HOME/code/petsc-3.3/
   export PETSC_ARCH=linux64


PETSc4py
########

::

   cd petsc/
   hg clone https://code.google.com/p/petsc4py 
   cd petsc4py
   patch -p1 < petsc4py.patch
   python setup.py build
   python setup.py install --user

More detailed instruction available in
http://code.google.com/p/petsc4py/source/browse/docs/source/install.rst

mpi4py
######

Again, if you can't use pip::

   curl -O https://mpi4py.googlecode.com/files/mpi4py-1.3.tar.gz
   tar -zxf mpi4py-1.3.tar.gz
   cd mpi4py-1.3.tar.gz
   python setup.py build
   python setup.py install --user

Detailed instructions in
http://mpi4py.scipy.org/docs/usrman/install.html

NumPy
#####

If you choose to build your own numpy, do that before building
petsc4py.

Thats's it! You can now follow the instructions to run the example.
