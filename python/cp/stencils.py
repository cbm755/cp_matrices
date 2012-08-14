"""
Spatial differentiation operator stencils for Closest Point Method
codes
"""
import numpy as np

def Laplacian_2nd(dim):
    """
    Second-order approximation to the Laplacian (the 5 pt stencil in
    2D)
    """
    if dim == 2:
        def DiffWeightsFcn(dx):
            """Weights for second-order approx to Laplacian"""
            if np.isscalar(dx):
                dx2 = dx*dx
                return np.array([-4., 1., 1., 1., 1.]) / dx2
            else:
                dx2 = dx[0]*dx[0]
                dy2 = dx[1]*dx[1]
                return np.array([-2./dx2-2./dy2, 1./dx2, 1./dx2,
                                 1./dy2, 1./dy2])
        # the points in the stencil
        DiffStencil = np.array([[ 0,  0],
                                [ 1,  0],
                                [-1,  0],
                                [ 0,  1],
                                [ 0, -1]])
        # stencil is typically a hypercross: longest arm of it is used
        # in calculating e.g., bandwidths
        DiffLongestArm = 1
    elif dim == 3:
        def DiffWeightsFcn(dx):
            if np.isscalar(dx):
                dx2 = dx*dx
                return np.array([-6., 1., 1., 1., 1., 1., 1.]) / dx2
            else:
                dx2 = dx[0]*dx[0]
                dy2 = dx[1]*dx[1]
                dz2 = dx[2]*dx[2]
                return np.array([-2./dx2-2./dy2-2./dz2, 1./dx2, 1./dx2,
                                 1./dy2, 1./dy2, 1./dz2, 1./dz2])
        DiffStencil = np.array([[ 0,  0,  0],
                                [ 1,  0,  0],
                                [-1,  0,  0],
                                [ 0,  1,  0],
                                [ 0, -1,  0],
                                [ 0,  0,  1],
                                [ 0,  0, -1]])
        DiffLongestArm = 1
    else:
        raise NotImplementedError('Dim {0} not implemented'.format(dim))
    return DiffWeightsFcn, DiffStencil, DiffLongestArm


def Laplacian_4th(dim):
    """
    Forth-order approximation to the Laplacian
    """
    if dim == 2:
        def DiffWeightsFcn(dx):
            d2 = dx*dx
            #return [ -5.0/d2, \
            #        (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
            #        (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2 ]
            print "** WARNING: hardcoded f96 stuff **"
            # TODO: this won't work if dx is integer
            #fl = type(dx)
            from numpy import float96
            fl = float96
            # TODO This doesn't exist in my Intel laptop
            # f1o12: float of 1/12
            f1o12 = fl(1)/fl(12)
            f4o3 = fl(4)/fl(3)
            # TODO should return an array
            return [-fl(5)/d2,
                     -f1o12/d2, f4o3/d2, f4o3/d2, -f1o12/d2,
                     -f1o12/d2, f4o3/d2, f4o3/d2, -f1o12/d2]
        DiffStencil = np.array([[ 0,   0],
                                [-2,   0],
                                [-1,   0],
                                [ 1,   0],
                                [ 2,   0],
                                [ 0,  -2],
                                [ 0,  -1],
                                [ 0,   1],
                                [ 0,   2]])
        DiffLongestArm = 2
    elif dim == 3:
        print 'WARNING: not tested yet since copy-paste from matlab code (but now cleaned up!)'
        def DiffWeightsFcn(dx):
            d2 = dx*dx
            return np.array([-15.0/2.0, -1.0/12.0, 4.0/3.0, 4.0/3.0,
                             -1.0/12.0, -1.0/12.0, 4.0/3.0, 4.0/3.0,
                             -1.0/12.0, -1.0/12.0, 4.0/3.0, 4.0/3.0,
                             -1.0/12.0]) / d2
        DiffStencil = np.array([[ 0,  0,  0],
                                [-2,  0,  0],
                                [-1,  0,  0],
                                [ 1,  0,  0],
                                [ 2,  0,  0],
                                [ 0, -2,  0],
                                [ 0, -1,  0],
                                [ 0,  1,  0],
                                [ 0,  2,  0],
                                [ 0,  0, -2],
                                [ 0,  0, -1],
                                [ 0,  0,  1],
                                [ 0,  0,  2]])
        DiffLongestArm = 2
    else:
        raise NotImplementedError('Dim {0} not implemented'.format(dim))
    return DiffWeightsFcn, DiffStencil, DiffLongestArm



def Biharmonic_2nd(dim):
    """
    WARNING: PROBABLY YOU DON'T WANT THIS FOR THE CLOSEST POINT
    METHOD.  It will be inconsistent.  Instead square the Laplacian
    matrix like in [Macdonald&Ruuth 2009].
    """
    if dim==2:
        def DiffWeightsFcn(dx):
            d4 = dx**4
            return np.array([-20., -1., 8., 8., -1., -1., 8., 8., -1.,
                             -2., -2., -2., -2.]) / d4
        DiffStencil = np.array([[ 0,  0],
                                [-2,  0],
                                [-1,  0],
                                [ 1,  0],
                                [ 2,  0],
                                [ 0, -2],
                                [ 0, -1],
                                [ 0,  1],
                                [ 0,  2],
                                [ 1,  1],
                                [ 1, -1],
                                [-1,  1],
                                [-1, -1]])
        DiffLongestArm = 2
    else:
        raise NotImplementedError('Dim {0} not implemented'.format(dim))
    return DiffWeightsFcn, DiffStencil, DiffLongestArm
