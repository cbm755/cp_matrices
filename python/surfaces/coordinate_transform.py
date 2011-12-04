''' coordinate transforms, cart2pol etc '''

from numpy import cos,sin,arctan2,hypot

def pol2cart(th, r, z=None):
    ''' Transform polar coordinates (th,r) into 2D cartesian coordinates.

    You can pass in z and it just returns it back.  TODO: should it
    make a copy?  How do deal with scalar case?

    Note: th, r in backwards "matlab order"

    Returns
    -------
    x
    y
    optionally z (a copy)
    '''

    if (z==None):
        return (r*cos(th), r*sin(th))
    else:
        return (r*cos(th), r*sin(th), z.copy())


def cart2pol(x, y, z=None):
    ''' Transform 2D cartesian coordinates into polar coordinates.

    Returns
    -------
    theta
    r
    optionally z (a copy)
    '''

    if (z==None):
        return (arctan2(y, x), hypot(x, y))
    else:
        return (arctan2(y, x), hypot(x, y), z.copy())


def cart2sph(x, y, z):
    ''' Stub: TODO implement '''
    raise NameError('not implemented yet')


def sph2cart(th,phi,r):
    ''' Stub: TODO implement '''
    # matlab uses: azimuth TH, elevation PHI, radius R
    raise NameError('not implemented yet')
