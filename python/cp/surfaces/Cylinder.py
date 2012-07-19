'''
    CP function for a cylinder in 3D
    TODO: generalize to add height's
'''

from Surface import ShapeWithBdy

class Cylinder(ShapeWithBdy):
    def __init__(self, radius=1.0, top=1.0, bottom=-1.0):
        self.radius = radius
        #self.top = top
        #self.bottom = bottom
        raise NameError('WIP')

    def closestPointToCartesian(self, xx):
        """TODO: top bottom not implemented"""
        x,y,z = xx
        th,r = cart2pol(x,y)
        (cpx,cpy) = pol2cart(th, self.radius)
        cp = xx.copy()
        cp[0] = cpx
        cp[1] = cpy
        cp[2] = z
        dist = norm(cp - xx, 2)
        return cp,dist,0,{}
