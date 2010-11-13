"""
Spatial differentiation operator stencils for Closest Point Method
codes
"""
from numpy import array as a

def Laplacian_2nd(dim):
    """
    Second-order approximation to the Laplacian (the 5 pt stencil in
    2D)
    """
    if (dim == 2):
        def f(dx):
            dx2 = dx*dx
            return [-4/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2]
        # a function of dx for the weights
        DiffWeightsFcn = f
        # the points in the stencil
        DiffStencil = [ a([ 0,  0]), \
                        a([ 1,  0]), \
                        a([-1,  0]), \
                        a([ 0,  1]), \
                        a([ 0, -1]) ]
        # stencil is typically a hypercross: longest arm of it is used
        # in calculating e.g., bandwidths
        DiffLongestArm = 1
    elif (dim == 3):
        def f(dx):
            dx2 = dx*dx
            return [-6/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0,  0]), \
                        a([ 1,  0,  0]), \
                        a([-1,  0,  0]), \
                        a([ 0,  1,  0]), \
                        a([ 0, -1,  0]), \
                        a([ 0,  0,  1]), \
                        a([ 0,  0, -1]) ]
        DiffLongestArm = 1
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)


def Laplacian_4th(dim):
    """
    Forth-order approximation to the Laplacian
    """
    if (dim == 2):
        def f(dx):
            d2 = dx*dx
            return [ -5.0/d2, \
                    (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
                    (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,   0]), \
                        a([-2,   0]), \
                        a([-1,   0]), \
                        a([ 1,   0]), \
                        a([ 2,   0]), \
                        a([ 0,  -2]), \
                        a([ 0,  -1]), \
                        a([ 0,   1]), \
                        a([ 0,   2]) ]
        DiffLongestArm = 2
    elif (dim == 3):
        print 'WARNING: not tested yet since copy-paste from matlab code'
        def f(dx):
            d2 = dx*dx
            return [ (-15.0/2.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0,  0]), \
		        a([-2,  0,  0]), \
		        a([-1,  0,  0]), \
		        a([ 1,  0,  0]), \
		        a([ 2,  0,  0]), \
		        a([ 0, -2,  0]), \
		        a([ 0, -1,  0]), \
		        a([ 0,  1,  0]), \
		        a([ 0,  2,  0]), \
		        a([ 0,  0, -2]), \
		        a([ 0,  0, -1]), \
		        a([ 0,  0,  1]), \
		        a([ 0,  0,  2]) ]
        DiffLongestArm = 2
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)



def Biharmonic_2nd(dim):
    """
    WARNING: PROBABLY YOU DON'T WANT THIS FOR THE CLOSEST POINT
    METHOD.  It is inconsistent.  Instead square the Laplacian matrix
    like in [Macdonald&Ruuth 2009].
    """
    if dim==2:
        def f(dx):
            d4 = dx*dx*dx*dx
            return [ -20/d4,  -1/d4,  8/d4,  8/d4,  -1/d4,  -1/d4,  8/d4,  8/d4,  -1/d4, \
                     -2/d4,  -2/d4,  -2/d4,  -2/d4 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0]), \
		        a([-2,  0]), \
		        a([-1,  0]), \
		        a([ 1,  0]), \
		        a([ 2,  0]), \
		        a([ 0, -2]), \
		        a([ 0, -1]), \
		        a([ 0,  1]), \
		        a([ 0,  2]), \
                        a([ 1,  1]), \
                        a([ 1, -1]), \
                        a([-1,  1]), \
                        a([-1, -1]) ]
        DiffLongestArm = 2
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)
