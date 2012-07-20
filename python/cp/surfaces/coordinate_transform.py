'''Coordinate transforms: cart2pol etc '''

from numpy import cos, sin, arctan2, hypot, sqrt

def pol2cart(th, r, z=None):
    '''Transform polar coordinates (th,r) into 2D cartesian coordinates.

    You can pass in z and it just returns it back.  TODO: should it
    make a copy?  How do deal with scalar case?

    Note: th, r in backwards "matlab order"

    Returns
    -------
    x
    y
    optionally z (a copy)
    '''

    if z is None:
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

    if z is None:
        return (arctan2(y, x), hypot(x, y))
    else:
        return (arctan2(y, x), hypot(x, y), z.copy())


def cart2sph(x, y, z):
    '''Transform cartesian coordinates to spherical.

    Like the matlab function! Returns (azimuth, elevation, radius).'''
    
    return (arctan2(y, x),
            arctan2(z, hypot(x, y)),
            sqrt(x**2 + y**2 + z**2))


def sph2cart(th, phi, r):
    '''Transform spherical coordinates to cartesian.

    Like the matlab function! (th is azimuth, phi is elevation)'''

    return (r * cos(phi) * cos(th),
            r * cos(phi) * sin(th),
            r * sin(phi))
