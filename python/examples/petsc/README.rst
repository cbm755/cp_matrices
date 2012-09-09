Colin's PETSc tests
-------------------

Serial matrix generate with PETSc solve
#######################################

First run heat_circle_serial.py::

    ipython
    run -i examples/petsc/heat_circle_serial.py

This will generate sparse matrices in SciPy, solve the problem using
those in serial.  It will also convert the scipy.sparse csr matrices to
PETSc Mat's and save them to disk.  There should now be several .dat
files and a .pickle file.

Then run heat_circle_petsc_solve.py::

    mpiexec -n 4 python examples/petsc/heat_circle_petsc_solve.py

or on my machine::

    PYTHONPATH="$PYTHONPATH:." mpiexec -n 4 python examples/petsc/heat_circle_petsc_solve.py

(this could be run on another machine, provided you copy the .pickle,
.dat and .dat.info files)


Output
######

The code will output the max norm of the error between the serial
solution and the PETSc solution.  It will indicate PASS/FAIL.


TODO
####

PETSc code with n processors seems to be at least n times as slow as
the serial code.  Perhaps its because I've only tested with explicit
matrix multiplies and relatively small matrices.

