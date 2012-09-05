"""
CPBar function.  Wraps an object and does the CPBar projection to it
if the closest point hit a boundary (bdy).

See section 5 in [Macdonald, Brandman, Ruuth] for an explanation of
\bar{cp}
"""
from Surface import ShapeWithBdy


class CPBar(ShapeWithBdy):
    def __init__(self, wrapped_shape):
        self.wrapped_shape = wrapped_shape

    def __getattr__(self, attr):
        """Get attributes not (re)implemented in this class from the
        wrapped shape.

        That is, if an attribute/method is implemented here, that's
        the one that'll be used."""
        return self.wrapped_shape.__getattribute__(attr)

    def closest_point(self, x):
        """ note: returns distance to the original cp """
        cpf = self.wrapped_shape.closest_point
        cpx, dist, bdy, others = cpf(x)
        if bdy==1 or bdy==2:
            y = x + 2*(cpx - x)
            cpx2, dist2, bdy2, others2 = cpf(y)
            others['origcp'] = cpx
            others['others2'] = others2
            if bdy2 != 0:
                #raise NameError('cpbar hit bdy!  What to do?')
                print 'cpbar hit bdy!  dist=', dist
                print x
                print dist, cpx, bdy
                print dist2, cpx2, bdy2
        else:
            cpx2 = cpx
        return cpx2, dist, bdy, others
