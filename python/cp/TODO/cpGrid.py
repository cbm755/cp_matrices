# make print a function as in python 3
#from __future__ import print_function

import numpy as np
from cpOps import findGridInterpBasePt, buildInterpWeights


# TODO: rename cpbaseptI should be cpbaseptIndex
class CPNode:
    """
    One grid point of a closest point tree data structure
    """
    def __init__(self, lsc, dx, level, basicPt):
        self.level = level
        self.gridpt = lsc
        self.dx = dx
        self.basicPt = basicPt
        self.children = []


    def __str__(self):
        return 'gridpt ' + str(self.gridpt) + \
            ' (level ' + str(self.level) + ')'


    def computeCP(self, cpfun):
        """
        Compute the CP for this node

        TODO: deal with "other" and bdy better
        """
        (cpx, dist, bdy, other) = cpfun(self.gridpt)
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
        # base 0, gridindex is 0 at relpt. needs to be consistent with
        # findGridInterpBasePt(), maybe this should be a function too
        self.gridIndex = ((self.gridpt - relpt) / self.dx).round().astype(int)
        # if (level >= 2):
        #     # sanity check, don't check on first level though
        #     # dx/10 is overkill: dx/2 maybe?  (must be larger than mach eps)
        eps = np.finfo(float).eps
        if (np.abs(self.gridIndex * self.dx + relpt - self.gridpt)).max() >= 10*eps:
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
                n = CPNode(lsc, dx, level+1, self.basicPt)
                n.computeCP(self.root.CPfun)
                n.computeGridIndex(self.root.basicPt)
                #maxdist = max(n.dist)
                maxdist = n.dist
                # TODO: this code should go elsewhere
                p = PAR_EXTSTENP
                dim = PAR_DIM
                arm = PAR_DiffLongestArm  # width of one "arm" of evolution stencil
                # Formula found on Ruuth & Merriman (for 2nd order centered
                # difference Laplacian, where arm=1). Beware, the grid spacing
                # is accounted in the comparison with maxdist, so in fact this
                # is \frac{\lambda}{\Delta x}
                lam = sqrt((dim-1.0)*((p+1)/2.0)**2 + (arm + (p+1)/2.0)**2)
                # Safety factor
                lam = 1.0001*lam
                if maxdist <= lam*dx:
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


    def plot(self, level=1, stop_at_level=None):
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

        if (level < stop_at_level):
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
        from time import time
        import stencils
        global PAR_DIM, PAR_EXTSTENP, PAR_EXTSTENWIDTH, PAR_EXTSTENSZ, \
            PAR_DiffLongestArm

        print "forming tree grid"
        st = time()
        #findGridInterpBasePt = findGridInterpBasePt
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
        PAR_EXTSTENSZ = PAR_EXTSTENWIDTH**dim

        if dim == 2:
            self.Dirs = [np.array([0,0]), np.array([1,0]), 
                         np.array([1,1]), np.array([0,1])]
        elif dim == 3:
            self.Dirs = [np.array([0,0,0]), np.array([1,0,0]), np.array([1,1,0]),
                         np.array([0,1,0]), np.array([0,1,1]), np.array([1,1,1]),
                         np.array([1,0,1]), np.array([0,0,1])]
        else:
            raise NotImplementedError('Dim {0} not implemented'.format(dim))

        # TODO: better way to default to Laplacian?  import stencils
        # outside the function?
        # TODO: want this more abstract, should just need to know the
        # maximum stencil here, might want to construct multiple operators
        # over a single grid.
        if diffInfo is None:
            diffInfo = stencils.Laplacian_2nd
        (self.DiffWeights, self.DiffStencil, PAR_DiffLongestArm) = diffInfo(dim)

        self.makeInterpStencil()

        self.Tree = []
        #self.Lists = [[]]*self.Levels  # no, same list many times
        self.Lists = [[] for i in xrange(self.Levels)]
        self.Grids = [dict() for i in xrange(self.Levels)]

        self.Lextend = [None for i in xrange(self.Levels)]
        self.Levolve = self.Lextend[:]
        self.Lghost = self.Lextend[:]
        self.Lnone = self.Lextend[:]

        self.populate()
        print "finished the tree grid, time = " + str(time() - st)

    def populate(self):
        n = CPNode(self.basicPt, self._dx, 1, self.basicPt)
        n.computeCP(self.CPfun)
        n.computeGridIndex(self.basicPt)
        self.addChild(n)
        n.subdivide(1)


    def makeInterpStencil(self):
        global PAR_EXTSTENWIDTH
        # TODO: sort out hte global dim with others
        dim = self.Dim
        self.InterpStencil = []
        # TODO: centralize
        #EXTSTENP = 3
        #EXTSTENWIDTH = EXTSTENP + 1
        if dim == 2:
            for i in range(PAR_EXTSTENWIDTH):
                for j in range (PAR_EXTSTENWIDTH):
                    self.InterpStencil.append(np.array([i,j]))
        elif dim == 3:
            for i in range(PAR_EXTSTENWIDTH):
                for j in range(PAR_EXTSTENWIDTH):
                    for k in range(PAR_EXTSTENWIDTH):
                        self.InterpStencil.append(np.array([i,j,k]))
        else:
            raise NameError('dim ' + str(dim) + ' not implemented')


    def plot(self, stop_at_level=5):
        for n in self.Tree:
            n.plot(stop_at_level=stop_at_level)



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

        print "Finding stencil sets at level " + str(level)
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
            except KeyError:
                print basept,relpt,basept*n.dx+relpt
                print n.dx
                print n,n.cp
                G[tuple(basept)]
                pass  #?


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
                        pass  #?
                    for off2 in self.DiffStencil:
                        gridind2 = basept + offset + off2
                        G[tuple(gridind2)].isExtendPt = True
        print "Found stencil sets in time = " + str(time() - st)


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
        print "Finding stencil sets at level", level
        st = time()
        Lextend = []
        Levolve = []
        Lghost = []
        Lnone = []
        G = self.Grids[level]
        for n in G.values():
            if n.isEvolvePt:
                if not n.isExtendPt:
                    raise ValueError('sanity fail')
                ind = len(Levolve)
                # n.IndexInLevolve = ind   # don't really need this
                # not a bug, its really Lextend here, see below
                n.IndexInLextend = ind
                Levolve.append(n)
            elif n.isExtendPt:
                if n.isEvolvePt:
                    raise ValueError('sanity fail')
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
        print "Found stencil sets in time =", time() - st


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
                n.interpweights = -buildInterpWeights(Xgrid, n.cp, n.dx, PAR_EXTSTENWIDTH)
            elif res == 'dirichlet_1st_order':
                # TODO: why bother with temp?  Just numpy.zeros?
                temp = buildInterpWeights(Xgrid, n.cp, n.dx, PAR_EXTSTENWIDTH)
                n.interpweights = [0]*len(temp)
            else:
                # includes "neumann_1st_order" and "neumann_2nd_order"
                n.interpweights = buildInterpWeights(Xgrid, n.cp, n.dx, PAR_EXTSTENWIDTH)
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



class CPFlatGrid(CPGrid):
    """
    A Closest Point grid

    todo: is inheritance best?
    """

    def __init__(self, name, cpfun, dim, lsc, dx, interp_degree=3, diffInfo=None):
        from time import time
        import stencils
        global PAR_DIM, PAR_EXTSTENP, PAR_EXTSTENWIDTH, PAR_EXTSTENSZ, \
            PAR_DiffLongestArm

        print "init grid"
        st = time()
        #findGridInterpBasePt = findGridInterpBasePt
        self.Name = name
        self.basicPt = lsc
        self.CPfun = cpfun
        self._dx = dx
        #self.Levels = levels
        levels = 1
        self.Dim = dim
        # this will be float64 or float96: make sure code accesses it
        # as neceesary
        self.floatType = type(dx)
        # TODO: do somthing special if dx is an int??

        PAR_DIM = dim
        PAR_EXTSTENP = interp_degree
        PAR_EXTSTENWIDTH = PAR_EXTSTENP+1
        PAR_EXTSTENSZ = PAR_EXTSTENWIDTH**dim

        if (dim == 2):
            self.Dirs = [np.array([0,0]), np.array([1,0]),
                         np.array([1,1]), np.array([0,1])]
        elif (dim == 3):
            self.Dirs = [np.array([0,0,0]), np.array([1,0,0]), np.array([1,1,0]),
                         np.array([0,1,0]), np.array([0,1,1]), np.array([1,1,1]),
                         np.array([1,0,1]), np.array([0,0,1])]
        else:
            raise NameError('dim ' + str(dim) + ' not implemented')

        # TODO: better way to default to Laplacian?  import stencils
        # outside the function?
        # TODO: want this more abstract, should just need to know the
        # maximum stencil here, might want to construct multiple operators
        # over a single grid.
        if diffInfo == None:
            diffInfo = stencils.Laplacian_2nd
        (self.DiffWeights, self.DiffStencil, PAR_DiffLongestArm) = diffInfo(dim)

        self.makeInterpStencil()

        self.Tree = [];
        #self.Lists = [[]]*self.Levels   # no, same list many times
        self.Lists = []
        for i in range(0, 1):
            self.Lists.append([])
        self.Grids = []
        for i in range(0, 1):
            self.Grids.append(dict())

        self.Lextend = []
        self.Levolve = []
        self.Lghost = []
        self.Lnone = []
        for i in range(0, 1):
            self.Lextend.append(None)
            self.Levolve.append(None)
            self.Lghost.append(None)
            self.Lnone.append(None)

        #self.populate()
        print "time = " + str(time() - st)


    def loadFromFile(self, fname):
        """
        Load the output from the fast C implementation of Ruuth's
        algorithm.
        """
        fp = np.float96
        f = open(fname, 'r')
        dx = self.dx
        relpt = self.center

        self.Grids = dict()

        class Node():
            pass

        #self.CPdd = dict()
        #self.CP = dict()
        #self.X = dict()

        for line in f:
            #print line,
            t = line.split()
            i = int(t[0])
            j = int(t[1])
            k = int(t[2])
            dd = fp(t[3])
            #cp = a([ fp(t[4]), fp(t[5]), fp(t[6]) ]);
            cp = fp(np.array([ t[4], t[5], t[6] ]));
            # TODO: issues here about loading float96's from file, see
            # my post to numpy maillist
            #x = fp(a([ t[7], t[8], t[9]]));
            x = np.array([ fp(t[7]), fp(t[8]), fp(t[9]) ]);
            x2 = np.array([i*dx+relpt[0], j*dx+relpt[1], k*dx+relpt[2]])
            #print "***", i, j, k, dd, cp, x-x2
            #self.CPdd[(i,j,k)] = dd
            #self.CP[(i,j,k)] = cp
            #self.X[(i,j,k)] = x2
            n = Node()
            n.dist = dd
            n.cp = cp
            n.gridpt = x2
            # TODO: better way to deal with this?
            n.whichBdy = 0
            self.Grids[(i,j,k)] = n
        # TODO: remove later, just for debuggin
        self.x = x
        self.x2 = x2
        f.close



if __name__ == "__main__":
    print "running as a script, import it instead"
