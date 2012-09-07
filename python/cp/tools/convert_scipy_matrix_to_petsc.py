A#def convert_scipy_matrix_to_petsc(A):

if (1==1):
    #import sys, petsc4py
    #petsc4py.init(sys.argv)
    from petsc4py import PETSc

    from scipy.sparse import find

    # TODO: OO way to break A into its list component form?
    Ai, Aj, As = find(A)

    #B = PETSc.Mat().create(PETSc.COMM_WORLD)
    B = PETSc.Mat().createAIJ((A.shape))
    #B.setSizes([m*n, m*n])
    #B.setSizes((A.shape))
    B.setFromOptions()
    Istart, Iend = B.getOwnershipRange()
    print (Istart, Iend)

    #A.setValue(0, 0, 3) # Insert a single value into matrix.
    #A.setValues([0, 1], [2, 3], [1, 1, 1, 1]) # Insert a 2x2 block of values into the matrix.

    # No!, but I'm sure there is a form like this
    #B.setValues(Ai, Aj, As)

    for i in xrange(0, A.shape[0]):
        I = Ai[i]
        J = Aj[i]
        val = As[i]
        #print (i,I,J,val)
        # TODO: this is silly, but I'll probably run in serial anyway...
        if I >= Istart and I < Iend:
        #B.setValue(I,J,val)
            B[I,J] = val

    B.assemblyBegin()
    B.assemblyEnd()

    x, y = B.getVecs()
    x.set(1)
    B.mult(x,y)

    viewer = PETSc.Viewer().createBinary('mymatrix.dat', 'w')
    viewer(B)
    #viewer = PETSc.Viewer().createBinary('vector-x.dat', 'w')
    #viewer(x)
    #viewer = PETSc.Viewer().createBinary('vector-y.dat', 'w')
    #viewer(y)

    # load
    viewer = PETSc.Viewer().createBinary('mymatrix.dat', 'r')
    B2 = PETSc.Mat().load(viewer)
    #viewer = PETSc.Viewer().createBinary('vector-x.dat', 'r')
    # = PETSc.Vec().load(viewer)
    #viewer = PETSc.Viewer().createBinary('vector-y.dat', 'r')
    #v = PETSc.Vec().load(viewer)

    # check
    assert B2.equal(B)
    #assert x.equal(u)
    #assert y.equal(v)

