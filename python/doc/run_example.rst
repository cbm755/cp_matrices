Run the examples
================

Based on some unknown issues, I can't plot the data by mayavi
on the cyclops.maths.ox.ac.uk. The solution now is to save
the solution into a file, copy the file to local machine and
plot the solution on the local machine.

Clone the repo
################

::

    #run on both of computing serve and local machine
    cd
    git clone https://github.com/cbm755/cp_matrices.git
    cd cp_matrices
    git checkout -t -b structure origin/structure
    git checkout structure

Prepare the mesh
################

::

    #run on both of computing server and local machine
    cd python
    cp /pathto/yourmesh.ply eight.ply

The file name should be eight.ply if you don't want to change the
code.  You can also put it in::

    python/cp/tests/data/eight.ply

(All of this is in a state of flux: no more ply files in the git repo
for now please.)

All above should be done in both of the computing server and the local
machine. The "Get the solution" part is of course time-consuming, and
only needs to done on the computing server.


Get the solution
###############

::

    #run on the computing server
    cd python
    mpiexec -n 24 python2.7 examples/cpm_heat_surface_ex.py
    # "24" should be the number of the processors.

After the long of computation, two files "gv.dat" and "gv.dat.info"
will be generated. Copy the gv.dat to local machine to plot the
solution.

Colin says: on my Fedora laptop, I have to do::

PYTHONPATH="$PYTHONPATH:." mpiexec -n 4 python examples/cpm_heat_surface_ex.py 

I (Jorge) think another way to solve it is by installing the cp
package. Something like::

    python setup.py install --user --record installed_files.txt

and then::

    cat installed_files.txt | xargs rm -rf

to remove the package if needed.

Copy the solution
###############

::

    #run on the local machine
    cd code/cp_matrices/python
    scp usrname@yourserve.edu:pathto/gv.dat .

Plot the solution
################

::

    # run on the local machine
    cd code/cp_matrices/python
    python2.7 examples/plot_mesh.py

I've successfully arrived here.

Notes and TODO
################

Colin: on my 32-bit machine I get an error about casting::

    Traceback (most recent call last):
      File "examples/cpm_heat_surface_ex.py", line 73, in <module>
        la,lv,gv,wv = band.createGLVectors()
     File "/home/cbm/work/cp_matrices-structure/python/cp/petsc/band.py", line 207, in createGLVectors
        self.SelectBlock()
      File "/home/cbm/work/cp_matrices-structure/python/cp/petsc/band.py", line 111, in SelectBlock
        self.ni2pi = PETSc.AO().createMapping(self.gindBlockWBand.getArray().astype(np.int64))
      File "AO.pyx", line 80, in petsc4py.PETSc.AO.createMapping (src/petsc4py.PETSc.c:135898)
      File "arraynpy.pxi", line 117, in petsc4py.PETSc.iarray_i (src/petsc4py.PETSc.c:6626)
      File "arraynpy.pxi", line 110, in petsc4py.PETSc.iarray (src/petsc4py.PETSc.c:6534)
    TypeError: array cannot be safely cast to required type
