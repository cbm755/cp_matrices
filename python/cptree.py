# make print a function as in python 3
#from __future__ import print_function

# TODO: rename cpbaseptI should be cpbaseptIndex

class CPTree:
    """
    One grid point of a closest point tree data structure
    """
    def __init__(self, lsc, dx, level):
        self.level = level
        self.gridpt = lsc
        self.dx = dx
        self.children = []


    def __str__(self):
        return 'gridpt ' + str(self.gridpt) + \
            ' (level ' + str(self.level) + ')'


    def computeCP(self, cpfun):
        (cpx, dist, bdy) = cpfun(self.gridpt)
        self.cp = cpx
        self.dist = dist
        self.whichBdy = bdy

    def computeGridIndex(self, relpt):
        """
        Compute the hash for this gridpt based on relpt.  This will be
        used in searching for points.  TODO: maybe hash it not the
        right word, its sort of an index.  NOTE: base is 0, so doesn't
        match my older matlab codes.
        """
        from numpy import finfo
        # base 0, gridindex is 0 at relpt. needs to be consistent with
        # findGridInterpBasePt(), maybe this should be a function too
        self.gridIndex = ((self.gridpt - relpt) / self.dx).round().astype(int)
        # if (level >= 2):
        #     # sanity check, don't check on first level though
        #     # dx/10 is overkill: dx/2 maybe?  (must be larger than mach eps)
        eps = finfo(float).eps
        if (abs(self.gridIndex * self.dx + relpt - self.gridpt)).max() >= 10*eps:
            raise NameError('gridindex failure')


    def subdivide(self, level):
        from math import sqrt
        global PAR_EXTSTENP, PAR_DIM, PAR_DiffLongestArm
        if (self.level <> level):
            error('level panic')

        maxLevels = self.root.Levels

        if level <= (maxLevels-1):
            #root = self.root;
            #flat = root.Lists{level};
            base = self.gridpt
            dx = self.dx / 2.0

            # TODO: no need to recompute CP for some point (but be
            # careful about caching, cpbar, etc)
            #for i in range(0, 4):
            #    lsc = base + self.root.Dirs[i]*dx
            #print self
            #print base
            #print dx
            for dir in self.root.Dirs:
                lsc = base + dir*dx
                n = CPTree(lsc, dx, level+1)
                n.computeCP(self.root.CPfun)
                n.computeGridIndex(self.root.basicPt)
                #maxdist = max(n.dist)
                maxdist = n.dist
                # TODO: this code should go elsewhere
                p = PAR_EXTSTENP
                dim = PAR_DIM
                arm = PAR_DiffLongestArm  # width of one "arm" of evolution stencil
                lam = sqrt( (dim-1.0)*((p+1)/2.0)**2 + (arm + (p+1)/2.0)**2)
                #print lam
                # safety factor
                lam = 1.0001*lam
                if (maxdist <= lam*dx):
                    self.addChild(n)
                    n.subdivide(level+1)


    def addChild(self, newNode):
        # TODO: implement these sanity checks
        #if isempty(parentNode.root)
        #warning('parent node not in a grid tree');

        #if ~isempty(newNode.Parent)
        #warning('node was already in a tree');

        #disconnect(newNode);
        newNode.parent = self
        newNode.root = self.root
        self.children.append(newNode)
        self.root.addToFlat(newNode)


    def plot(self, level=1):
        import pylab

        styles = ['r-o', 'b-x', 'g-s', 'm-+', 'k-.']*3
        #randcol = rand(3,1);
        #%offset = node.Dx / 10 * rand(2,1);
        #offset = [0;0];

        #pylab.plot([self.gridpt[0]], [self.gridpt[1]], 'rx')

        #pylab.plot([self.gridpt[0],self.cp[0]], [self.gridpt[1],self.cp[1]], \
        #               styles[level-1])

        pylab.plot([self.gridpt[0]], [self.gridpt[1]], styles[level-1])
        #pylab.plot([self.cp[0]], [self.cp[1]], styles[level-1])

        for n in self.children:
            n.plot(level+1)


    def findNeighbours(self):
        print 'finding neighbours of ' + str(self)
        # first strategy is just search the flat grids (might be slow)
        level = self.level
        F = self.root.Lists[level-1]

        n = self
        print "self's gridindex:" + str(n.gridIndex)
        for dir in self.root.DiffStencil:
            #print dir
            nbGridIndex = n.gridIndex + dir
            #print nbGridIndex
            foundIt = False
            index = -42
            # TODO:
            # for i,n2 in enumerate(F):
            for i in range(0, len(F)):
                n2 = F[i]
                if (nbGridIndex == n2.gridIndex).all():
                    print nbGridIndex, n2.gridIndex
                    foundIt = True
                    index = i
            if (foundIt):
                print 'found it, index=' + str(index)
                nb = F[index]
            else:
                print 'can''t find'
                # TODO: assign n,s,e,w to None
        #print F

        # TEMP
        if (level <= 2):
            for n in self.children:
                n.findNeighbours()



class CPGrid:
    """
    A Closest Point grid
    
    TODO: remain Lists and Grids as _RawLists and _RawGrids, then make
    Lists and Grids the reasonable grids only.  TODO: how to find
    reasonable?  Its the first grid that contains the full bandwidth,
    but how do we know that?
    
    """


    def __init__(self, name, cpfun, dim, lsc, dx, interp_degree=3, diffInfo=None, levels=6):
        #import numpy
        #a = numpy.array
        from numpy import array as a
        from time import time
        import stencils
        global PAR_DIM, PAR_EXTSTENP, PAR_EXTSTENWIDTH, PAR_EXTSTENSZ, \
            PAR_DiffLongestArm

        print "forming tree grid"
        st = time()
        self.Name = name
        self.basicPt = lsc
        self.CPfun = cpfun
        self._dx = dx
        self.Levels = levels
        self.Dim = dim
        # this will be float64 or float96: make sure code accesses it
        # as neceesary
        self.floatType = type(dx)
        # TODO: do somthing special if dx is an int??

        PAR_DIM = dim
        PAR_EXTSTENP = interp_degree
        PAR_EXTSTENWIDTH = PAR_EXTSTENP+1
        PAR_EXTSTENSZ = PAR_EXTSTENWIDTH**PAR_DIM

        if (dim == 2):
            self.Dirs = [a([0,0]), a([1,0]), a([1,1]), a([0,1])]
        elif (dim == 3):
            self.Dirs = [a([0,0,0]), a([1,0,0]), a([1,1,0]), a([0,1,0]), \
                         a([0,1,1]), a([1,1,1]), a([1,0,1]), a([0,0,1])]
        else:
            raise NameError('dim ' + str(dim) + ' not implemented')

        # TODO: better way to default to Laplacian?  import stencils
        # outside the function?
        if diffInfo == None:
            diffInfo = stencils.Laplacian_2nd
        (self.DiffWeights, self.DiffStencil, PAR_DiffLongestArm) = diffInfo(dim)

        self.makeInterpStencil()

        self.Tree = [];
        #self.Lists = [[]]*self.Levels   # no, same list many times
        self.Lists = []
        for i in range(0, self.Levels):
            self.Lists.append([])
        self.Grids = []
        for i in range(0, self.Levels):
            self.Grids.append(dict())

        self.Lextend = []
        self.Levolve = []
        self.Lghost = []
        self.Lnone = []
        for i in range(0, self.Levels):
            self.Lextend.append(None)
            self.Levolve.append(None)
            self.Lghost.append(None)
            self.Lnone.append(None)

        self.populate()
        print "finished the tree grid, time = " + str(time() - st)

    def populate(self):
        n = CPTree(self.basicPt, self._dx, 1)
        n.computeCP(self.CPfun)
        n.computeGridIndex(self.basicPt)
        self.addChild(n);
        n.subdivide(1);


    def makeInterpStencil(self):
        from numpy import array as a
        global PAR_EXTSTENWIDTH
        # TODO: sort out hte global dim with others
        dim = self.Dim
        self.InterpStencil = []
        # TODO: centralize
        #EXTSTENP = 3
        #EXTSTENWIDTH = EXTSTENP + 1
        if dim == 2:
            for j in range(0, PAR_EXTSTENWIDTH):
                for i in range (0, PAR_EXTSTENWIDTH):
                    self.InterpStencil.append( a([i,j]) )
        elif dim == 3:
            for k in range(0, PAR_EXTSTENWIDTH):
                for j in range(0, PAR_EXTSTENWIDTH):
                    for i in range(0, PAR_EXTSTENWIDTH):
                        self.InterpStencil.append( a([i,j,k]) )
        else:
            raise NameError('dim ' + str(dim) + ' not implemented')


    def plot(self):
        for n in self.Tree:
            n.plot()



    def addChild(self, newNode):
        #TODO
        #if ~isempty(newNode.Parent)
        #  warning('node was already in a tree');

        # TODO?
        #newNode.disconnect();

        # TODO: is this how to do it?
        newNode.parent = self
        newNode.root = self
        self.Tree.append(newNode)
        self.addToFlat(newNode)


    def addToFlat(self, newNode):
        self.Lists[newNode.level-1].append(newNode)
        # TODO: could do this a posterior with dict.fromkeys()
        self.Grids[newNode.level-1][tuple(newNode.gridIndex)] = newNode


    def findNeighbours(self):
        for n in self.Tree:
            n.findNeighbours()


    def findStencilSetsFaster(self, level, extraPts=[]):
        """
        Find the stencil sets [Macdonald & Ruuth 2008].  Here we use a
        faster version: mark the CP base points in the grid.  Then
        just iterate over those in a separate loop
        """
        from time import time
        global PAR_EXTSTENP

        print "finding stencil sets at level " + str(level)
        st = time()
        #F = self.Lists[level]
        G = self.Grids[level]
        # find the the CP, then find the points in the extsten
        # surrouding it, mark each of them
        relpt = self.basicPt
        # first, make them all false
        for n in G.values():
            n.isEvolvePt = False
            n.isExtendPt = False
            n.isCPBasePt = False
        # essentially the algorithm in [MR2008], but here we mark each
        # CP Base Point and then iterate over those in a separate
        # loop.  I confirmed it gives the same thing as the slower
        # version and much faster
        for n in G.values():
            cp = n.cp
            basept = findGridInterpBasePt(cp,n.dx,relpt,PAR_EXTSTENP)
            n.cpbaseptI = basept
            try:
                G[tuple(basept)].isCPBasePt = True
            except (KeyError):
                print basept,relpt,basept*n.dx+relpt
                print n.dx
                print n,n.cp
                G[tuple(basept)]
                pass


        for n in G.values():
            if n.isCPBasePt == True:
                basept = n.gridIndex
                for offset in self.InterpStencil:
                    gridindex = basept + offset
                    try:
                        G[tuple(gridindex)].isEvolvePt = True
                    except (KeyError):
                        print basept,offset,gridindex
                        G[tuple(basept)]
                        pass
                    for off2 in self.DiffStencil:
                        gridind2 = basept + offset + off2
                        G[tuple(gridind2)].isExtendPt = True
        print "found stencil sets in time = " + str(time() - st)


    def findStencilSetsSlower(self, level, extraPts=[]):
        """Find the stencil sets [Macdonald & Ruuth 2008]"""
        from time import time
        global PAR_EXTSTENP

        print "finding stencil sets at level " + str(level)
        st = time()
        #F = self.Lists[level]
        G = self.Grids[level]
        # find the the CP, then find the points in the extsten
        # surrouding it, mark each of them
        relpt = self.basicPt
        # first, make them all false
        for n in G.values():
            n.isEvolvePt = False
            n.isExtendPt = False
        # now, follow the algorithm in [MR2008]
        for n in G.values():
            cp = n.cp
            basept = findGridInterpBasePt(cp,n.dx,relpt,PAR_EXTSTENP)
            n.cpbaseptI = basept
            for offset in self.InterpStencil:
                gridindex = basept + offset
                # for each extrapt, find pts in the extsten and mark them
                #try:
                G[tuple(gridindex)].isEvolvePt = True
                #except (KeyError):
                #    raise NameError('not enough points in grid')
                for off2 in self.DiffStencil:
                    gridind2 = basept + offset + off2
                    #    try:
                    G[tuple(gridind2)].isExtendPt = True
                    #    except (KeyError):
                    #        raise NameError('not enough points in grid')
        print "found stencil sets in time = " + str(time() - st)


    findStencilSets = findStencilSetsFaster


    def buildListsFromStencilSets(self, level):
        """ """
        from time import time
        print "finding stencil sets at level " + str(level)
        st = time()
        Lextend = []
        Levolve = []
        Lghost = []
        Lnone = []
        G = self.Grids[level]
        for n in G.values():
            if n.isEvolvePt:
                if not n.isExtendPt: raise NameError('sanity fail')
                ind = len(Levolve)
                # n.IndexInLevolve = ind   # don't really need this
                # not a bug, its really Lextend here, see below
                n.IndexInLextend = ind
                Levolve.append(n)
            elif n.isExtendPt:
                if n.isEvolvePt: raise NameError('sanity fail')
                Lghost.append(n)
            else:
                #print "  not in either list"
                Lnone.append(n)
        Lextend = list(Levolve)
        Lextend.extend(Lghost)
        # We can change a quadratic search time to linear in a later
        # part of the code by storing the indices into Lextend here.
        # We found the index for the ones in Levolve already, and the
        # index matches that of Lextend.  Now we just need to add the
        # indices for the ones in Lghost.
        lenLevolve = len(Levolve)
        for i,n in enumerate(Lghost):
            n.IndexInLextend = lenLevolve + i
        self.Lextend[level] = Lextend
        self.Levolve[level] = Levolve
        self.Lghost[level] = Lghost
        self.Lnone[level] = Lnone
        print "found stencil sets in time = " + str(time() - st)


    def bdyfcn_null(bdy):
        return 0

    def findStencilsOnStencilSets(self, level, bdyfcn=bdyfcn_null):
        from time import time
        print "icpm preparartion: diff"
        st = time()
        Lextend = self.Lextend[level]
        Levolve = self.Levolve[level]
        #Lghost = self.Lghost[level]
        G = self.Grids[level]
        for n in Levolve:
            n.diffpts = []
            gi = n.gridIndex
            for offsets in self.DiffStencil:
                gii = gi + offsets
                nn = G[tuple(gii)]
                #mm = Lextend.index(nn)
                # Previous works but is quadratic (this is significant
                # in practice).  Instead we find the indices just
                # after we build the stencils and cache them.
                mm = nn.IndexInLextend
                n.diffpts.append(mm)
        print "  time=" + str(time()-st)

        print "icpm preparartion: interp"
        st = time()
        for n in Lextend:
            n.interppts = []
            cpbaseptI = n.cpbaseptI
            Xgrid = G[tuple(cpbaseptI)].gridpt
            res = bdyfcn(n.whichBdy)
            if res == 'dirichlet_2nd_order':
                n.interpweights = -buildInterpWeights(Xgrid, n.cp, n.dx)
            elif res == 'dirichlet_0th_order':
                temp = buildInterpWeights(Xgrid, n.cp, n.dx)
                n.interpweights = [0]*len(temp)
            else:
                # includes "neumann_0th_order" and "neumann_2nd_order"
                n.interpweights = buildInterpWeights(Xgrid, n.cp, n.dx)
            #gi = n.gridIndex
            for offsets in self.InterpStencil:
                # TODO: not right: find basept, UPDATE: this is
                # basept, should be correct now
                gii = cpbaseptI + offsets
                nn = G[tuple(gii)]
                #mm = Levolve.index(nn)
                # Previous works but is quadratic, see above in the
                # diff code.  Levolve and Lextend have the same
                # indices so we use index in extend here
                mm = nn.IndexInLextend
                n.interppts.append(mm)
        print "  time=" + str(time()-st)


    def buildDiffMatrix(self, level):
        """
        generate the matrix D
        """
        from math import log10,ceil
        import numpy
        import scipy.sparse
        from time import time

        Levolve = self.Levolve[level]
        Lextend = self.Lextend[level]
        dx = Levolve[0].dx
        progout = 10 ** (ceil(log10(len(Levolve)))-1)
        print "building D"
        st=time()
        diffWeights = self.DiffWeights(dx)
        # TODO: float96 fixes?
        #D = scipy.sparse.lil_matrix((len(Levolve),len(Lextend)))
        #D = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
        #                                 dtype=type(diffWeights[0]) )
        D = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                         dtype=self.floatType )
        for i,n in enumerate(Levolve):
            if i % progout == 0:
                print "  D row " + str(i)
            for s,j in enumerate(n.diffpts):
                D[i,j] = diffWeights[s]
        print "  D row " + str(len(Levolve))
        print "elapsed time = " + str(time()-st)
        return D


    def buildDiffMatrixTempDxDyTest(self, level):
        """
        generate the matrices Dx-, Dy-
        Hardcoded for 2D
        """
        from math import log10,ceil
        import numpy
        import scipy.sparse
        from time import time

        if (self.Dim != 2):
            raise NameError('this routine hardcoded for dimension 2')
        Levolve = self.Levolve[level]
        Lextend = self.Lextend[level]
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
                                          dtype=self.floatType )
        Dxf = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                          dtype=self.floatType )
        Dyb = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                          dtype=self.floatType )
        Dyf = scipy.sparse.lil_matrix( (len(Levolve),len(Lextend)), \
                                          dtype=self.floatType )

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


    def buildExtensionMatrix(self, level):
        """
        generate the matrix E
        """
        from math import ceil,log10
        import scipy.sparse
        from time import time

        Levolve = self.Levolve[level]
        Lextend = self.Lextend[level]
        dx = Levolve[0].dx
        progout = 10 ** (ceil(log10(len(Lextend)))-1)
        print "building E"
        st=time()
        #TODO: floattype fixes
        #E = scipy.sparse.lil_matrix((len(Lextend),len(Levolve)))
        E = scipy.sparse.lil_matrix( (len(Lextend),len(Levolve)), \
                                         dtype=self.floatType )

        for i,n in enumerate(Lextend):
            if i % progout == 0:
                print "  E row " + str(i)
            for s,j in enumerate(n.interppts):
                E[i,j] = n.interpweights[s]
        print "  E row " + str(len(Lextend))
        print "elapsed time = " + str(time() - st)
        return E




def LagrangeWeights1D(xg, x, dx, N):
    """
    1D Lagrange interpolation weights
    xg is the base grid point
    """
    from numpy import zeros, array as a
    from scipy import comb
    # Is exactly on a grid point, then return binary weights
    # TODO: here they are returned as integers: maybe should be fp,
    # also need to find if its float64 or float96.
    # TODO: Nick Hale mentioned they just do the calculation and than
    # catch the NaN/Inf exception and return the binary weights; said
    # that was faster.
    w = zeros(N,dtype=type(dx))
    for j in range(0, N):
        if x == (xg+j*dx):
            w[j] = 1
            return w

    # TODO: comb very slow and floating point: better to cache the
    # values I think, see prun
    #for i in range(0,N):
    #    w[i] = (-1)**i * comb(N-1,i)
    # TODO: slightly faster to use some sort of switch statement here?
    if   (N == 1): ww = [1]
    elif (N == 2): ww = [1, -1]
    elif (N == 3): ww = [1, -2, 1]
    elif (N == 4): ww = [1, -3, 3, -1]
    elif (N == 5): ww = [1, -4, 6, -4, 1]
    elif (N == 6): ww = [1, -5, 10, -10, 5, -1]
    elif (N == 7): ww = [1, -6, 15, -20,  15,  -6,   1]
    elif (N == 8): ww = [1, -7, 21, -35,  35, -21,   7,  -1]
    else: raise NameError('need to hardcode more weights')

    # careful about types here, maybe double sth?
    for j in range(0,N):
        w[j] = ww[j] / ( x - (xg + j*dx) )
    w = w / sum(w);
    return w


def LagrangeWeights1DSlow(xg, x, dx, N):
    """
    1D Lagrange interpolation weights
    xg is the base grid point

    This version is slow: and maybe doesn't work with float96
    """
    from numpy import zeros, array as a
    from scipy import comb
    # Is exactly on a grid point, then return binary weights
    # TODO: here they are returned as integers: maybe should be fp,
    # also need to find if its float64 or float96
    w = zeros(N)
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


def buildInterpWeights(Xgrid, X, dx):
    """
    build the interpolation weights
    Xg is the B-pt (the base grid point below X)
    X is the point to be evaluated at
    """
    global PAR_DIM, PAR_EXTSTENSZ, PAR_EXTSTENWIDTH
    from numpy import finfo, zeros


    xweights = LagrangeWeights1D(Xgrid[0], X[0], dx, PAR_EXTSTENWIDTH)
    yweights = LagrangeWeights1D(Xgrid[1], X[1], dx, PAR_EXTSTENWIDTH)
    if (PAR_DIM == 3):
        zweights = LagrangeWeights1D(Xgrid[2], X[2], dx, PAR_EXTSTENWIDTH)

    s = 0
    sum1 = 0.0*dx
    # need a type here because of float96
    extWeights = zeros( PAR_EXTSTENSZ, dtype=type(dx) )
    if (PAR_DIM == 2):
        # loop in the same order as elsewhere and compute weights as
        # products of above
        for j in range(0,PAR_EXTSTENWIDTH):
            for i in range(0,PAR_EXTSTENWIDTH):
                extWeights[s] = xweights[i]*yweights[j]
                sum1 = sum1 + extWeights[s]
                s = s + 1
                if s > PAR_EXTSTENSZ:
                    raise NameError('too many points in stencil')
    elif (PAR_DIM == 3):
        for k in range(0,PAR_EXTSTENWIDTH):
            for j in range(0,PAR_EXTSTENWIDTH):
                for i in range(0,PAR_EXTSTENWIDTH):
                    extWeights[s] = xweights[i]*yweights[j]*zweights[k]
                    sum1 = sum1 + extWeights[s]
                    s = s + 1
                    if s > PAR_EXTSTENSZ:
                        raise NameError('too many points in stencil')
    else:
        raise NameError('that dimension not implemented')

    # check here depends on what type of float
    # TODO: is this sanity check expensive?
    eps = finfo(type(dx)).eps
    if abs(sum1 - 1.0) > 50*eps:
        print 'weight problem'
        print extWeights
        print sum1
        raise NameError('weight problem')

    return extWeights


def findGridInterpBasePt(x,dx,relpt,EXTSTENP):
    """find the basepoint for x"""
    # p=0: ==B==
    # p=1:   B====x
    # p=2:   B  ==x==  x
    # p=3:   B    x====x    x
    # p=4:   B    x  ==x==  x    x
    # p=5:   B    x    x====x    x    x
    # etc
    # This index is base 0 (i.e., B=(0,0,...) is the bottom left of
    # the whole grid)
    from numpy import floor,round
    if EXTSTENP % 2 == 0:  # even
        #I = round( (PAR.NPOINTS-1)/ (PAR.DomB-PAR.DomA) * (x-PAR.DomA)  ) + 1;
        I = round( (x-relpt) / dx ).astype(int) + 1
        B = I - (EXTSTENP)/2 - 1
    else:  # odd
        #I = floor( (PAR.NPOINTS-1) * (x-PAR.DomA) / (PAR.DomB-PAR.DomA) ) + 1;
        I = floor( (x-relpt) / dx ).astype(int) + 1
        B = I - (EXTSTENP - 1)/2 - 1
    return B


def LinearDiagonalSplitting(D, E):
    """
    Compute the DEstab matrix
    """
    import scipy.sparse

    usz, lsz = D.shape

    # old scipy (<7?) version
    if (1==0):
        #I = scipy.sparse.lil_eye((10,10),k=0)
        #Inonsq = scipy.sparse.lil_eye((usz,usz+gsz),k=0)
        Ddiagv = scipy.sparse.extract_diagonal(D)
        #Ddiagpadded = scipy.sparse.lil_diags([Ddiagv], [0], (usz, usz+gsz))
        #Ddiagpadded = scipy.sparse.spdiags(Ddiagv, 0, usz, usz+gsz)
        #Ddiagpadded = scipy.sparse.lil_matrix(shape=(usz,usz+gsz))
        Ddiag = scipy.sparse.lil_eye((usz,usz))
        Ddiag.setdiag(Ddiagv)
        Ddiagpadded = scipy.sparse.lil_eye((usz,usz+gsz))
        Ddiagpadded.setdiag(Ddiagv)
        DE = Ddiag + (D - Ddiagpadded)*E

    Ddiagv = D.diagonal()
    #Ddiag    = scipy.sparse.lil_diags([Ddiagv], [0], (usz,usz))
    #Ddiagpad = scipy.sparse.lil_diags([Ddiagv], [0], (usz,usz+gsz))
    Ddiag    = scipy.sparse.spdiags(Ddiagv, 0, usz,usz)
    Ddiagpad = scipy.sparse.spdiags(Ddiagv, 0, usz,lsz)
    DE = Ddiag + (D - Ddiagpad)*E
    return DE

DEstab = LinearDiagonalSplitting


def buildEPlotMatrix(g, level, Vertices, VertexBpt = None):
    # uband, usz, gband, gsz, Vertices, VertexBpt, PAR):
    """
    Interpolate to get values at Vertices
    """
    import numpy
    import scipy.sparse
    from time import time
    from math import log10,ceil
    global PAR_EXTSTENP

    # TODO: VertexBpt could cache the Bpts for each vertex

    # TODO: we could find two matrices: E and some more rows, then
    # take the product with tihs?  this is to deal with the problem
    # that these points in Vertices may need grid

    st = time()
    print "building Eplot"
    Lextend = g.Lextend[level]
    Levolve = g.Levolve[level]
    G = g.Grids[level]

    dx = Levolve[0].dx
    relpt = g.basicPt
    N = Vertices.shape[0]
    progout = 10 ** (ceil(log10(N))-1)

    EPlot = scipy.sparse.lil_matrix( (N,len(Lextend)), dtype=type(dx) )

    for i in range(0,N):
        if i % progout == 0:
            print "  Eplot row " + str(i)

        #   % find floor of x,y,z in uband and the CP, these are the only two
        #   % things we need to construct this row of E
        x = Vertices[i,:]
        # TODO: could cache the baseptI
        #Bpt = VertexBpt[i]
        #if uband.isBpt[Bpt] != 1:
        #    raise NameError("should be a Bpt!!")
        xbaseptIndex = findGridInterpBasePt(x, dx, relpt, PAR_EXTSTENP)

        Xgrid = G[tuple(xbaseptIndex)].gridpt
        interpWeights = buildInterpWeights(Xgrid, x, dx)

        for s,offsets in enumerate(g.InterpStencil):
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

    return EPlot



if __name__ == "__main__":
    print "running as a script, import it instead"
