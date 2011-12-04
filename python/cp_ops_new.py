def buildDiffMatrix(g, Levolve, Lextend):
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
    diffWeights = Levolve[0].root.DiffWeights(dx)
    stencilsize = len(diffWeights)
    #D = lil_matrix( (len(Levolve),len(Lextend)), dtype=type(dx) )
    # make empty lists for i,j and a_{ij}
    ii = [];  jj = [];  aij = []
    #for i,n in enumerate(Levolve):
    for i in range(0, len(g.band)):
        if i % progout == 0:  print "  D row " + str(i)
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
