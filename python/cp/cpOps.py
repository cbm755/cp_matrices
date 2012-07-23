"""
closest point matrix operators

TODO: these access too many internals of the CPGrid obj, just pass
Levol and Lext?  Still working on the separation of code here...

TODO: split out barycentric Lagrange interp
"""

import numpy as np

def buildDiffMatrixNoDX(Levolve, Lextend):
    """
    generate the matrix D
    """
    from math import log10,ceil
    from scipy.sparse import coo_matrix
    from time import time

    progout = 10 ** (ceil(log10(len(Levolve)))-1)
    print "building D"
    st=time()

    dx = Levolve[0].dx
    diffWeights = Levolve[0].root.DiffWeights(1)
    stencilsize = len(diffWeights)
    #D = lil_matrix( (len(Levolve),len(Lextend)), dtype=type(dx) )
    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    for i,n in enumerate(Levolve):
        if i % progout == 0:
            print "  D row " + str(i)
        ## Slow approach using lil_matrix
        #for s,j in enumerate(n.diffpts):
        #    D[i,j] = diffWeights[s]
        ## Fast coo_matrix approach
        ii.extend([i]*stencilsize)
        jj.extend(n.diffpts)
        aij.extend(diffWeights)
    D = coo_matrix( (aij,(ii,jj)), shape=(len(Levolve),len(Lextend)), dtype=type(dx) )
    print "  D row " + str(len(Levolve))
    print "elapsed time = " + str(time()-st)
    return D.tocsr()

def buildDiffMatrix(Levolve, Lextend):
    """
    Generate the matrix D.
    """
    from math import log10, ceil
    from scipy.sparse import coo_matrix
    from time import time

    progout = 10 ** (ceil(log10(len(Levolve)))-1)
    print "building D"
    st=time()

    dx = Levolve[0].dx
    diffWeights = Levolve[0].root.DiffWeights(dx)
    stencilsize = len(diffWeights)
    #D = lil_matrix( (len(Levolve),len(Lextend)), dtype=type(dx) )
    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    for i, n in enumerate(Levolve):
        if i % progout == 0:  print "  D row " + str(i)
        ## Slow approach using lil_matrix
        #for s,j in enumerate(n.diffpts):
        #    D[i,j] = diffWeights[s]
        ## Fast coo_matrix approach
        ii.extend([i]*stencilsize)
        jj.extend(n.diffpts)
        aij.extend(diffWeights)
    D = coo_matrix((aij,(ii,jj)), shape=(len(Levolve),len(Lextend)),
                   dtype=type(dx))
    print "D row " + str(len(Levolve))
    print "elapsed time =", time()-st
    return D.tocsr()


def _buildDiffMatrix_lil_slow(Levolve, Lextend):
    """
    Generate the matrix D using lil_matrix
    
    Deprecated, slower than coo_matrix
    """
    from math import log10,ceil
    from scipy.sparse import lil_matrix
    from time import time
    progout = 10 ** (ceil(log10(len(Levolve)))-1)
    print "building D"
    st=time()
    dx = Levolve[0].dx
    diffWeights = Levolve[0].root.DiffWeights(dx)
    #diffWeights = Grid.DiffWeights(dx)
    D = lil_matrix( (len(Levolve),len(Lextend)), dtype=type(dx) )
    for i,n in enumerate(Levolve):
        if i % progout == 0:
            print "  D row " + str(i)
        for s,j in enumerate(n.diffpts):
            D[i,j] = diffWeights[s]
    print "  D row " + str(len(Levolve))
    print "elapsed time = " + str(time()-st)
    return D.tocsr()



def _buildDiffMatrix_depreciated(Grid, level):
    """
    Generate the matrix D, old version that needs to know Grid
    internals
    """
    from math import log10,ceil
    import numpy
    from scipy.sparse import lil_matrix
    from time import time

    Levolve = Grid.Levolve[level]
    Lextend = Grid.Lextend[level]
    dx = Levolve[0].dx
    progout = 10 ** (ceil(log10(len(Levolve)))-1)
    print "building D"
    st=time()
    diffWeights = Grid.DiffWeights(dx)
    # TODO: float96 fixes?
    #D = lil_matrix( (len(Levolve),len(Lextend)) )
    #D = lil_matrix( (len(Levolve),len(Lextend)), dtype=type(diffWeights[0]) )
    D = lil_matrix( (len(Levolve),len(Lextend)), dtype=Grid.floatType )
    for i,n in enumerate(Levolve):
        if i % progout == 0:
            print "  D row " + str(i)
        for s,j in enumerate(n.diffpts):
            D[i,j] = diffWeights[s]
    print "  D row " + str(len(Levolve))
    print "elapsed time = " + str(time()-st)
    return D.tocsr()


def buildDiffMatrixTempDxDyTest(Grid, level):
    """
    generate the matrices Dx-, Dy-
    Hardcoded for 2D
    """
    from math import log10,ceil
    import scipy.sparse
    from time import time

    if (Grid.Dim != 2):
        raise NameError('this routine hardcoded for dimension 2')
    Levolve = Grid.Levolve[level]
    Lextend = Grid.Lextend[level]
    dx = Levolve[0].dx
    progout = 10 ** (ceil(log10(len(Levolve)))-1)
    print "building D"
    st=time()
    # the points in the stencil
    # DiffStencil = [ a([ 0,  0]), \
    #                 a([ 1,  0]), \
    #                 a([-1,  0]), \
    #                 a([ 0,  1]), \
    #                 a([ 0, -1]) ]
    diffWxb = [ 1.0/dx, 0.0,   -1.0/dx, 0.0,    0.0]
    diffWxf = [-1.0/dx, 1.0/dx, 0.0,    0.0,    0.0]
    diffWyb = [ 1.0/dx, 0.0,    0.0,    0.0,   -1.0/dx]
    diffWyf = [-1.0/dx, 0.0,    0.0,    1.0/dx, 0.0]


    Dxb = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                      dtype=Grid.floatType )
    Dxf = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                      dtype=Grid.floatType )
    Dyb = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                      dtype=Grid.floatType )
    Dyf = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                      dtype=Grid.floatType )

    for i,n in enumerate(Levolve):
        if i % progout == 0:
            print "  Dx,Dy row " + str(i)
        for s,j in enumerate(n.diffpts):
            Dxb[i,j] = diffWxb[s]
            Dxf[i,j] = diffWxf[s]
            Dyb[i,j] = diffWyb[s]
            Dyf[i,j] = diffWyf[s]
    print "  D row " + str(len(Levolve))
    print "elapsed time = " + str(time()-st)
    return (Dxb,Dxf,Dyb,Dyf)



def buildExtensionMatrix(Levolve, Lextend):
    r"""
    Generate the matrix E

    Notes:
    In the diff op case, the weights are fixed, here they depend on
    which node we're at.
    """
    from math import ceil, log10
    from scipy.sparse import coo_matrix
    from time import time

    dx = Levolve[0].dx
    stencilsize = len(Levolve[0].interppts)
    progout = 10 ** (ceil(log10(len(Lextend)))-1)
    print "building E"
    st=time()

    #E = lil_matrix( (len(Lextend),len(Levolve)), dtype=type(dx) )
    # make empty lists for i,j and e_{ij}
    ii = [];  jj = [];  aij = []
    for i,n in enumerate(Lextend):
        if i % progout == 0:  print "  E row " + str(i)
        ## Slowest: using lil_matrix
        #for s,j in enumerate(n.interppts):
        #    E[i,j] = n.interpweights[s]
        ## Fast: using coo_matrix
        #for s,j in enumerate(n.interppts):
        #    ii.append(i)
        #    jj.append(j)
        #    dd.append(n.interpweights[s])
        ## Fastest: vectorized and using coo_matrix (roughly twice as fast)
        ii.extend([i]*stencilsize)
        jj.extend(n.interppts)
        aij.extend(n.interpweights)
    #TODO: does this work with float96?
    E = coo_matrix((aij,(ii,jj)), shape=(len(Lextend),len(Levolve)),
                   dtype=type(dx))
    print "  E row " + str(len(Lextend))
    print "elapsed time = " + str(time() - st)
    return E.tocsr()


def _buildExtensionMatrix_lil_slow(Levolve, Lextend):
    r"""
    Generate the matrix E

    Notes:
    In the diff op case, the weights are fixed, here they depend on
    which node we're at.

    TODO: would it be faster to collect $i$, $j$ and $a_{ij}$ then
    build all at once?
    """
    from math import ceil,log10
    from scipy.sparse import lil_matrix
    from time import time

    #Levolve = Grid.Levolve[level]
    #Lextend = Grid.Lextend[level]
    dx = Levolve[0].dx
    progout = 10 ** (ceil(log10(len(Lextend)))-1)
    print "building E"
    st=time()
    #TODO: floattype fixes
    E = lil_matrix( (len(Lextend),len(Levolve)), dtype=type(dx) )

    for i,n in enumerate(Lextend):
        if i % progout == 0:
            print "  E row " + str(i)
        for s,j in enumerate(n.interppts):
            E[i,j] = n.interpweights[s]
    print "  E row " + str(len(Lextend))
    print "elapsed time = " + str(time() - st)
    return E.tocsr()


def _buildExtensionMatrix_depreciated(Grid, level):
    """
    generate the matrix E, old version that uses Grid directly
    """
    from math import ceil,log10
    from scipy.sparse import lil_matrix
    from time import time

    Levolve = Grid.Levolve[level]
    Lextend = Grid.Lextend[level]
    dx = Levolve[0].dx
    progout = 10 ** (ceil(log10(len(Lextend)))-1)
    print "building E"
    st=time()
    #TODO: floattype fixes
    #E = lil_matrix( (len(Lextend),len(Levolve)) )
    E = lil_matrix( (len(Lextend),len(Levolve)), dtype=Grid.floatType )

    for i,n in enumerate(Lextend):
        if i % progout == 0:
            print "  E row " + str(i)
        for s,j in enumerate(n.interppts):
            E[i,j] = n.interpweights[s]
    print "  E row " + str(len(Lextend))
    print "elapsed time = " + str(time() - st)
    return E.tocsr()



def findGridInterpBasePt(x,dx,relpt,p):
    """
    Find the "base grid point" for a point x.  This is best explained
    in the diagram below.

    x: is the interpolation point, must lie inside the === signs below
    p: degree interpolation (N-1 point interp)
    relpt: a reference point, corresponding to (0,0) in your grid (can be a vector)
    dx: grid spacing (can be a vector)
    B: index to the "basepoint", the lower left corner of an interpolation
stencil hypercube.  It is a 2D/3D/etc index, measured relative to
'relpt'

    p=0: ==B==
    p=1:   B====x
    p=2:   B  ==x==  x
    p=3:   B    x====x    x
    p=4:   B    x  ==x==  x    x
    p=5:   B    x    x====x    x    x
    etc
    This index is base 0 (i.e., B=(0,0,...) is the bottom left of the
    whole grid)

    TODO: processes one point: should accept an array for x (maybe it already works?)
    """
    if p % 2 == 0:  # even
        I = np.round( (x-relpt) / dx ).astype(int) + 1
        B = I - p/2 - 1
    else:  # odd
        I = np.floor( (x-relpt) / dx ).astype(int) + 1
        B = I - (p - 1)/2 - 1
    return B



def LagrangeWeights1D(xg, x, dx, N):
    """1D Lagrange interpolation weights for equidistant points.

    Based on Berrut & Trefethen "Barycentric Lagrange Interpolation".

    Input
    -----
    xg : base grid point
    x : array of lenght M, or scalar
        Points for which the interpolation weights will be calculated.
    dx : grid spacing
    N : number of interpolation points (interpolation degree + 1)

    Output
    ------
    w : array of shape (M,N) if x is an array, else array of shape (N,)
        Barycentric weights.

    Usage
    -----
    Let f be a function such that f(xg + i \cdot dx) = f_i for 0 \leq
    i \leq N-1

    To interpolate in point x, call
    np.dot(LagrangeWeights1D(xg, x, dx, N), [f_i for i in range(N)])
    """
    #from scipy import comb
    # Is exactly on a grid point, then return binary weights
    # TODO: here they are returned as integers: maybe should be fp,
    # also need to find if its float64 or float96.

    # To accept an array x, we have to catch NaNs (else, it needs
    # quite some bookkeeping)
    
    # TODO: comb very slow and floating point: better to cache the
    # values I think, see prun
    # These are the weights w_j for equidistant nodes
    # (5.1) in Berrut & Trefethen "Barycentric Lagrange Interpolation"
    #for i in range(0,N):
    #    w[i] = (-1)**i * comb(N-1,i)
    try:
        ww = {1:np.array([1]),
              2:np.array([1, -1]),
              3:np.array([1, -2, 1]),
              4:np.array([1, -3, 3, -1]),
              5:np.array([1, -4, 6, -4, 1]),
              6:np.array([1, -5, 10, -10, 5, -1]),
              7:np.array([1, -6, 15, -20,  15,  -6,   1]),
              8:np.array([1, -7, 21, -35,  35, -21,   7,  -1]),
              }[N]
    except KeyError:
        raise ValueError('need to hardcode more weights')

    scalar = np.isscalar(x)
    x = np.atleast_1d(x)
    # Maybe from __future__ import division just to be sure?
    # Barycentric formula (4.2) in Berrut & Trefethen
    # To avoid getting RuntimeWarning
    np.seterr(invalid='ignore', divide='ignore')
    w = ww / (x[:, np.newaxis] - (xg + np.arange(N) * dx))
    w /= np.sum(w, axis=1)[:, np.newaxis]
    np.seterr(invalid='warn', divide='warn')  # back to default settings
    w[np.isnan(w)] = 1
    # 15% faster using bottleneck here:
    # import bottleneck as bn
    # bn.replace(w, np.nan, 1.)
    if scalar:
        return w[0]
    else:
        return w


def LagrangeWeights1DSlow(xg, x, dx, N):
    """
    1D Lagrange interpolation weights
    xg is the base grid point

    This version is slow: and maybe doesn't work with float96
    """
    from scipy.misc import comb
    # Is exactly on a grid point, then return binary weights
    # TODO: here they are returned as integers: maybe should be fp,
    # also need to find if its float64 or float96
    w = np.zeros(N)
    for j in range(0, N):
        if x == (xg+j*dx):
            w[j] = 1
            return w

    for i in range(0,N):
        w[i] = (-1)**i * comb(N-1,i)

    # careful about types here, maybe double sth?
    for j in range(0,N):
        w[j] = w[j] / ( x - (xg + j*dx) )
    w = w / sum(w);
    return w


def buildInterpWeights(Xgrid, X, dx, EXTSTENWIDTH):
    """
    build the interpolation weights
    Xgrid is the B-pt (the base grid point below X)
    X is the point to be evaluated at

    TODO: clean up EXTSTENWIDTH and/or document (its degree+1)
    """
    dim = len(X)

    EXTSTENSZ = EXTSTENWIDTH**dim

    if np.isscalar(dx):
        dxv = [dx]*3
    else:
        dxv = dx

    xweights = LagrangeWeights1D(Xgrid[0], X[0], dxv[0], EXTSTENWIDTH)
    yweights = LagrangeWeights1D(Xgrid[1], X[1], dxv[1], EXTSTENWIDTH)
    if (dim == 3):
        zweights = LagrangeWeights1D(Xgrid[2], X[2], dxv[2], EXTSTENWIDTH)

    #print extWeights.dtype, xweights.dtype, yweights.dtype
    if dim == 2:
        # loop in the same order as elsewhere and compute weights as
        # products of above
        extWeights = (yweights[:, np.newaxis] * xweights[np.newaxis, :]).ravel()
    elif dim == 3:
        extWeights = (zweights[:, np.newaxis, np.newaxis] *  # z varies the slowest, so put it in the first dimension
                      yweights[np.newaxis, :, np.newaxis] *
                      xweights[np.newaxis, np.newaxis, :]).ravel()  # x varies the fastest, put it in the last dimension
    else:
        raise NotImplementedError('Dimension not implemented yet.')

    sum1 = np.sum(extWeights, dtype=type(dxv[0]))
    # check here depends on what type of float
    # TODO: is this sanity check expensive?
    # Yes it is: calculating sum1 and checking this takes about the same time
    # as calculating extWeights given {x, y, z}weights
    eps = np.finfo(type(dxv[0])).eps
    if abs(sum1 - 1.0) > 50*eps:
        print extWeights
        print sum1
        raise ValueError('Weight problem')

    return extWeights


def LinearDiagonalSplitting(D, E):
    """
    Compute the DEstab matrix.

    Implements stable modification of the implicit Closest Point
    Mehtod procedure, see formula (2.8) in [ICPM].
    """
    import scipy.sparse

    usz, lsz = D.shape
    Ddiagv = D.diagonal()
    Ddiag = scipy.sparse.spdiags(Ddiagv, 0, usz, usz)
    Ddiagpad = scipy.sparse.spdiags(Ddiagv, 0, usz, lsz)
    # sparse matrices, so this is matrix-matrix multiplication, not elementwise.
    DE = Ddiag + (D - Ddiagpad)*E
    return DE

# Alias, matches language in [Macdonald & Ruuth 2009]
DEstab = LinearDiagonalSplitting


#def buildEPlotMatrix(g, level, Vertices, interp_degree, VertexBpt = None):
def buildEPlotMatrix(G, Levolve, Lextend, Points, interp_degree, PointsBpt = None):
    """
    Interpolate to get values at Points

    TODO: currently, the degree here myst match the "main" one, but
    this could be changed (see e.g., g.InterpStencil())
    """
    from scipy.sparse import coo_matrix
    from time import time
    from math import log10, ceil

    if interp_degree == 0:
        return buildEPlotMatrixNN(G, Levolve, Lextend, Points)
    # TODO
    #elif interp_degree != Levolve[0].root.interp_degree:
    #    raise NameError()


    relpt = Levolve[0].basicPt
    # TODO: this is duplicated info... or maybe not if this degree can
    # differ
    EXTSTENP = interp_degree
    EXTSTENWIDTH = EXTSTENP+1
    interpStencil = Levolve[0].root.InterpStencil

    # TODO: VertexBpt could cache the Bpts for each vertex

    # TODO: we could find two matrices: E and some more rows, then
    # take the product with tihs?  this is to deal with the problem
    # that these points in Vertices may need grid

    st = time()
    print "building Eplot"
    #Lextend = g.Lextend[level]
    #Levolve = g.Levolve[level]
    #G = g.Grids[level]

    dx = Levolve[0].dx
    #relpt = g.basicPt
    N = Points.shape[0]
    progout = 10 ** (ceil(log10(N))-1)

    #EPlot = lil_matrix( (N,len(Lextend)), dtype=type(dx) )

    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    for i in range(0,N):
        if i % progout == 0:  print "  Eplot row " + str(i)
        #   % find floor of x,y,z in uband and the CP, these are the only two
        #   % things we need to construct this row of E
        x = Points[i,:]
        # TODO: could cache the baseptI
        #Bpt = VertexBpt[i]
        #if uband.isBpt[Bpt] != 1:
        #    raise NameError("should be a Bpt!!")
        xbaseptIndex = findGridInterpBasePt(x, dx, relpt, EXTSTENP)

        Xgrid = G[tuple(xbaseptIndex)].gridpt
        interpWeights = buildInterpWeights(Xgrid, x, dx, EXTSTENWIDTH)

        for s,offsets in enumerate(interpStencil):
                gii = xbaseptIndex + offsets
                nn = G[tuple(gii)]
                #mm = Levolve.index(nn)
                # Previous works but is quadratic, see above in the
                # diff code.  Levolve and Lextend have the same
                # indices so we use index in extend here
                mm = nn.IndexInLextend
                # TOOD: how about we just make it Eplot*E*u?
                # TODO: should we check if we find it in evolve here?
                #if (mm >= len(Levolve)):
                #    print "about to fail on row " + str(i) + ", mm=" \
                #        + str(mm) + ", len(Levolve)=" + str(len(Levolve))
                # TODO: may not even be in extend!?  can this even happen?
                if (mm >= len(Lextend)):
                    raise NameError("Outside Lextend, report")

                #EPlot[i,mm] = interpWeights[s]
                # ii.extend([i]*stencilsize)
                # jj.extend(n.diffpts)
                # aij.extend(diffWeights)
                ii.append(i)
                jj.append(mm)
                aij.append(interpWeights[s])
    EPlot = coo_matrix( (aij,(ii,jj)), shape=(N,len(Lextend)), dtype=type(dx) )
    print "elapsed time = " + str(time()-st)

    # LIL/COO matrices are slow in most use cases, convert to CSR for
    # faster matrix-vector products (something like 24x faster on
    # 10000x10000 in my tests)
    return EPlot.tocsr()


def buildEPlotMatrixNN(G, Levolve, Lextend, Points):
    """
    Interpolate to get values at Points using nearest neighbours (fast, sparse)
    """
    from scipy.sparse import coo_matrix
    from time import time
    from math import log10,ceil

    raise NameError('TODO')

    dim = len(Levolve[0].gridpt)
    relpt = Levolve[0].basicPt
    # TODO: this is duplicated info... or maybe not if this degree can
    # differ
    EXTSTENP = interp_degree
    EXTSTENWIDTH = EXTSTENP+1
    interpStencil = Levolve[0].root.InterpStencil

    # TODO: VertexBpt could cache the Bpts for each vertex

    # TODO: we could find two matrices: E and some more rows, then
    # take the product with tihs?  this is to deal with the problem
    # that these points in Vertices may need grid

    st = time()
    print "building Eplot"
    #Lextend = g.Lextend[level]
    #Levolve = g.Levolve[level]
    #G = g.Grids[level]

    dx = Levolve[0].dx
    #relpt = g.basicPt
    N = Points.shape[0]
    progout = 10 ** (ceil(log10(N))-1)

    #EPlot = lil_matrix( (N,len(Lextend)), dtype=type(dx) )

    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    for i in range(0,N):
        if i % progout == 0:  print "  Eplot row " + str(i)
        x = Points[i,:]
        xbaseptIndex = findGridInterpBasePt(x, dx, relpt, EXTSTENP)

        Xgrid = G[tuple(xbaseptIndex)].gridpt
        interpWeights = buildInterpWeights(Xgrid, x, dx, EXTSTENWIDTH)

        for s,offsets in enumerate(interpStencil):
                gii = xbaseptIndex + offsets
                nn = G[tuple(gii)]
                #mm = Levolve.index(nn)
                # Previous works but is quadratic, see above in the
                # diff code.  Levolve and Lextend have the same
                # indices so we use index in extend here
                mm = nn.IndexInLextend
                # TOOD: how about we just make it Eplot*E*u?
                # TODO: should we check if we find it in evolve here?
                #if (mm >= len(Levolve)):
                #    print "about to fail on row " + str(i) + ", mm=" \
                #        + str(mm) + ", len(Levolve)=" + str(len(Levolve))
                # TODO: may not even be in extend!?  can this even happen?
                if (mm >= len(Lextend)):
                    raise NameError("Outside Lextend, report")

                #EPlot[i,mm] = interpWeights[s]
                # ii.extend([i]*stencilsize)
                # jj.extend(n.diffpts)
                # aij.extend(diffWeights)
                ii.append(i)
                jj.append(mm)
                aij.append(interpWeights[s])
    EPlot = coo_matrix( (aij,(ii,jj)), shape=(N,len(Lextend)), dtype=type(dx) )
    print "elapsed time = " + str(time()-st)

    return EPlot.tocsr()



def buildEPlotMatrix_old(G, Levolve, Lextend, Points, interp_degree, PointsBpt = None):
    """
    Interpolate to get values at Points

    TODO: currently, the degree here myst match the "main" one, but
    this could be changed (see e.g., g.InterpStencil())

    TODO: passing in relpt seems like a bad idea.  G should know about this.
    """
    #import numpy
    from scipy.sparse import lil_matrix
    from time import time
    from math import log10,ceil

    #dim = g.Dim
    dim = len(Levolve[0].gridpt)
    relpt = Levolve[0].basicPt
    # TODO: this is duplicated info... or maybe not if this degree can
    # differ
    EXTSTENP = interp_degree
    EXTSTENWIDTH = EXTSTENP+1
    interpStencil = Levolve[0].root.InterpStencil

    # TODO: VertexBpt could cache the Bpts for each vertex

    # TODO: we could find two matrices: E and some more rows, then
    # take the product with tihs?  this is to deal with the problem
    # that these points in Vertices may need grid

    st = time()
    print "building Eplot"
    #Lextend = g.Lextend[level]
    #Levolve = g.Levolve[level]
    #G = g.Grids[level]

    dx = Levolve[0].dx
    #relpt = g.basicPt
    N = Points.shape[0]
    progout = 10 ** (ceil(log10(N))-1)

    EPlot = lil_matrix( (N,len(Lextend)), dtype=type(dx) )

    for i in range(0,N):
        if i % progout == 0:
            print "  Eplot row " + str(i)

        #   % find floor of x,y,z in uband and the CP, these are the only two
        #   % things we need to construct this row of E
        x = Points[i,:]
        # TODO: could cache the baseptI
        #Bpt = VertexBpt[i]
        #if uband.isBpt[Bpt] != 1:
        #    raise NameError("should be a Bpt!!")
        xbaseptIndex = findGridInterpBasePt(x, dx, relpt, EXTSTENP)

        Xgrid = G[tuple(xbaseptIndex)].gridpt
        interpWeights = buildInterpWeights(Xgrid, x, dx, EXTSTENWIDTH)

        for s,offsets in enumerate(interpStencil):
                gii = xbaseptIndex + offsets
                nn = G[tuple(gii)]
                #mm = Levolve.index(nn)
                # Previous works but is quadratic, see above in the
                # diff code.  Levolve and Lextend have the same
                # indices so we use index in extend here
                mm = nn.IndexInLextend
                # TOOD: how about we just make it Eplot*E*u?
                # TODO: should we check if we find it in evolve here?
                #if (mm >= len(Levolve)):
                #    print "about to fail on row " + str(i) + ", mm=" \
                #        + str(mm) + ", len(Levolve)=" + str(len(Levolve))
                # TODO: may not even be in extend!?  can this even happen?
                if (mm >= len(Lextend)):
                    raise NameError("Outside Lextend, report")

                EPlot[i,mm] = interpWeights[s]

    print "elapsed time = " + str(time()-st)

    # LIL matrices are slow in most use cases, convert to CSR for
    # faster matrix-vector products (something like 24x faster on
    # 10000x10000 in my tests)
    return EPlot.tocsr()
