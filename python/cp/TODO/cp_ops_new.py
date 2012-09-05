import stencils
from numpy import array as a
from cpOps import findGridInterpBasePt, buildInterpWeights

def buildDiffMatrix(g):
    """
    generate the matrix D
    TODO: right now defaults to Laplacian but this should be a general
    function
    """
    from math import log10,ceil
    from scipy.sparse import coo_matrix
    from time import time


    progout = 10 ** (ceil(log10(len(g.band)))-1)
    print "building D using loops"
    st=time()

    # get the stencil info
    (f,stencil,armlen) = stencils.Laplacian_2nd(dim=2)
    diffWeights = f(g.dx)
    stencilsize = len(diffWeights)

    (warn1,warn2,warn3,warn4) = (0,0,0,0)

    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    #for i,n in enumerate(Levolve):
    for c in range(0, len(g.band)):
        if c % progout == 0:  print "  D row " + str(c)
        ij = g.ij[c:c+1]
        ii.extend([c]*stencilsize)
        # no banding:
        #ii.extend([g.band[c]]*stencilsize)
        for s,pt in enumerate(stencil):
            ij2 = ij + pt
            # TODO, clean up
            if (ij2[0][0] < 0):
                warn1 = warn1 + 1
                ij2[0][0] = g.nx - 1
            if (ij2[0][1] < 0):
                warn2 = warn2 + 1
                ij2[0][1] = g.ny - 1
            if (ij2[0][0] >= g.nx):
                warn3 = warn3 + 1
                ij2[0][0] = 0
            if (ij2[0][1] >= g.ny):
                warn4 = warn4 + 1
                ij2[0][1] = 0
            c2 = g.sub2ind(ij2)
            #print (c,s,pt,ij,ij2,c2)
            jj.extend(c2)
            aij.extend([diffWeights[s]])

    #print (len(aij), len(ii), len(jj))
    #print (len(g.band), g.nx*g.ny)
    #print (max(aij), max(ii), max(jj))
    #print (min(aij), min(ii), min(jj))

    # TODO: support float96?  dtype=
    D = coo_matrix( (aij,(ii,jj)), shape=(len(g.band), g.nx*g.ny) )
    print "  D row " + str(len(g.band))
    if warn1 > 0:
        print "  warning: periodic BC applied on left %d times" % warn1
    if warn2 > 0:
        print "  warning: periodic BC applied on bottom %d times" % warn2
    if warn3 > 0:
        print "  warning: periodic BC applied on right %d times" % warn2
    if warn4 > 0:
        print "  warning: periodic BC applied on top %d times" % warn2
    print "  elapsed time = " + str(time()-st)
    return D.tocsr()


def buildDiffMatrixFast(g):
    """
    generate the matrix D
    TODO: right now defaults to Laplacian but this should be a general
    function
    """
    from math import log10,ceil
    from scipy.sparse import coo_matrix
    from time import time
    from numpy import array as a

    progout = 10 ** (ceil(log10(len(g.band)))-1)
    print "building D using vectors"
    st=time()

    # get the stencil info
    (f,stencil,armlen) = stencils.Laplacian_2nd(dim=2)
    diffWeights = f(g.dx)
    stencilsize = len(diffWeights)

    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []

    ij = g.ij
    #I = g.band
    #I2 = g.sub2ind(ij)

    BCwarn = 0

    nx = g.nx
    ny = g.ny
    for s,pt in enumerate(stencil):
        ij2 = ij + pt
        # TODO, clean up, do this BC stuff in a separate function
        # (ideally in a cleaner way, with some warning when this
        # happens)
        M1 = ij2 < a([0,-42])
        M2 = ij2 < a([-42,0])
        M3 = ij2 > a([nx-1,2*ny])
        M4 = ij2 > a([2*nx,ny-1])
        ij2 = (M1) * (nx-1) + (~M1) * ij2
        ij2 = (M2) * (ny-1) + (~M2) * ij2
        ij2 = (M3) * 0 + (~M3) * ij2
        ij2 = (M4) * 0 + (~M4) * ij2
        BCwarn = BCwarn | M1.any() | M2.any() | M3.any() | M4.any()

        I2 = g.sub2ind(ij2)
        #ii.extend(I)  # BUG!  this won't be banded in the rows
        ii.extend(range(0,len(g.band)))
        jj.extend(I2)
        aij.extend([diffWeights[s]]*len(g.band))

    #print (len(aij), len(ii), len(jj))
    #print (len(g.band), g.nx*g.ny)
    #print (max(aij), max(ii), max(jj))
    #print (min(aij), min(ii), min(jj))

    # TODO: support float96?  dtype=
    D = coo_matrix( (aij,(ii,jj)), shape=(len(g.band), g.nx*g.ny) )
    if BCwarn > 0:
        print "  warning: periodic BC applied at least once"
    print "  elapsed time = " + str(time()-st)
    return D.tocsr()



def buildExtensionMatrix(g, xy, degreep=3):
    r"""
    Generate the matrix E

    Notes:
    In the diff op case, the weights are fixed, here they depend on
    which node we're at.
    """
    from math import ceil,log10
    from scipy.sparse import coo_matrix
    from time import time
    from numpy import zeros, isscalar

    # TODO: dim hardcoded to 2 in some places in this function
    dim = xy.shape[1]
    dx = g.dx
    if (isscalar(dx)):
        mytype = type(dx)
    else:
        mytype = type(dx[0])

    relpt = a((g.x1d[0], g.y1d[0]))
    #stencilsize = len(Levolve[0].interppts)
    stencilsize = (degreep+1)**dim
    progout = 10 ** (ceil(log10(len(g.band)))-1)
    print "building E"
    st=time()

    # make empty lists the sparse matrix
    # TODO: is it faster to use numpy arrays here?
    ii = [];  jj = [];  aij = []

    # how many points to interpolate
    N = xy.shape[0]

    # TODO: order here must match code elsewhere: this is a very bad idea
    stencil = zeros( (stencilsize,dim) , dtype=int)
    c = 0
    for j in range(0,degreep+1):
        for i in range(0,degreep+1):
            stencil[c][0] = i
            stencil[c][1] = j
            c = c+1

    # can do all of them at once
    #bpij = findGridInterpBasePt(xy, dx, relpt, degreep)
    for c in range(0, N):
        if c % progout == 0:  print "  E row " + str(c)
        x = xy[c]
        bpij = findGridInterpBasePt(x, dx, relpt, degreep)
        bpx = relpt + bpij * dx
        interpweights = buildInterpWeights(bpx, x, dx, degreep+1)
        interppts = g.sub2ind(bpij + stencil)

        ii.extend([c]*stencilsize)
        jj.extend(interppts)
        aij.extend(interpweights)
    #TODO: does this work with float96?
    E = coo_matrix( (aij,(ii,jj)), shape=(len(g.band),g.nx*g.ny), dtype=mytype )

    # TODO: return this info as well:
    from numpy import unique
    band2 = unique(jj)
    print (len(jj),len(band2),len(g.band))
    print "  E row " + str(len(g.band))
    print "elapsed time = " + str(time() - st)
    #return E.tocsr()
    #E2 = E.tocsr()
    #EE = E2[:,band2]
    return (E.tocsr(), band2)
