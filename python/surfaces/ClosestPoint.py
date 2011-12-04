#import numpy
from numpy import array as a

class ClosestPoint():
    def __init__(self):
        print 'base class constructor'

    def closestPointToCartesian(self, x):
        raise NameError('should be implemented in subclass')
        #cp = numpy.zeros(x.shape)
        #dist = numpy.linalg.norm(x-cp,2)
        #return cp,dist,bdy,{}

    def cpwrap(self, x):
        return self.closestPointToCartesian(self, x)
    #def closestPointToPolar(self, x):


    def getBB(self):
        """
        Return a bounding box for the object.
        Could be overridden by subclasses.
        """
        return self._bb


    def viztest(self, extra_bb=1.3):
        if self._dim == 2:
            self._viztest2d(extra_bb=extra_bb)
        elif self._dim == 3:
            self._viztest3d(extra_bb=extra_bb)
        else:
            raise NameError('Only 2d and 3d viz tests implemented')

    def _viztest2d(self, extra_bb):
        import pylab
        import numpy
        plot = pylab.plot

        if (self._hasParam):
            X,Y = self.ParamGrid()
            plot(X,Y,'k-')

        a,b = self.getBB()
        ll = (b+a)/2 - (extra_bb)*(b-a)/2
        rr = (b+a)/2 + (extra_bb)*(b-a)/2
        #xx = numpy.random.uniform(bb[0], bb[1], (100,2))
        for i in range(0,100):
            # TODO: here we assume the object lives in [-2,2]^2
            # TODO: maybe each object could record a boundingbox
            #x = 4*numpy.random.random((2)) - 2
            x = numpy.random.uniform( ll, rr )
            cp,dist,bdy,other = self.cp(x)
            col = [0.4, 0.4, 0.4]
            plot([x[0], cp[0]], [x[1],cp[1]], '-', color=col)
            plot([cp[0]], [cp[1]], 'o', color=col)
        pylab.show()

    def _viztest3d(self, extra_bb):
        """
        3D vizualization of CPRep for obj
        """
        from enthought.mayavi import mlab
        import numpy as np
        f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(640,640))
        mlab.clf()

        if self._hasParam:
            L = self.ParamGrid()
            x,y,z = L
            #print x,y,z
            #s = mlab.mesh(xb, yb, zb, scalars=real((Eplot_bulb*E*uxx).reshape(xb.shape)))
            # TODO: mayavi bug, opacity<1
            s = mlab.mesh(x, y, z, scalars=z, opacity=1.0)

        a,b = self.getBB()
        ll = (b+a)/2 - (extra_bb)*(b-a)/2
        rr = (b+a)/2 + (extra_bb)*(b-a)/2
        for i in range(0,40):
            x = np.random.uniform( ll, rr )
            cp,dist,bdy,other = self.cp(x)
            colt = (0.5,0.5,0.5)
            op = 0.3
            l = mlab.plot3d([x[0],cp[0]], [x[1],cp[1]], [x[2], cp[2]], color=colt, opacity=op, tube_radius=0.1)
        #mlab.title('3d viz test')
        mlab.show()




class ShapeWithBdy(ClosestPoint):
    def _viztest2d(self, extra_bb):
        import pylab
        import numpy
        plot = pylab.plot

        if (self._hasParam):
            X,Y = self.ParamGrid()
            plot(X,Y,'k-')

        a,b = self.getBB()
        ll = (b+a)/2 - (extra_bb)*(b-a)/2
        rr = (b+a)/2 + (extra_bb)*(b-a)/2
        for i in range(0,100):
            x = numpy.random.uniform( ll, rr )
            cp,dist,bdy,other = self.cp(x)
            if bdy==0:
                col = [0.4, 0.4, 0.4]
            elif bdy==1:
                col = 'r'
            elif bdy==2:
                col = 'b'
            else:
                # TODO:
                print bdy
                raise NameError('should do something for other bdy values')
            plot([x[0], cp[0]], [x[1],cp[1]], '-', color=col)
            plot([cp[0]], [cp[1]], 'o', color=col)
        axis('scaled')
        pylab.show()


    def _viztest3d(self):
        """
        3D vizualization of CPRep for obj with boundary
        """
        from enthought.mayavi import mlab
        import numpy as np
        f = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(500,700))
        mlab.clf()

        #s = mlab.mesh(xb, yb, zb, scalars=real(evec_b.reshape(xb.shape)))
        #l = mlab.plot3d(xs, ys, zs, real(evec_s))
        #mlab.title(str(ii) + ' ew=' + str(eval), size=0.2)

        #mlab.show()
        #mlab.savefig('b_horn' + str(ii) + '_' + str(eval) + '.png')
        #s = mlab.mesh(xb, yb, zb, scalars=real((Eplot_bulb*E*uxx).reshape(xb.shape)))
        #l = mlab.plot3d(xs, ys, zs, real((Eplot_stem*E*uxx).reshape(xs.shape)))

        #(x1,y1),(x2,y2),(x3,y3) = mesh2d(resolution=3)

        # TODO: is has param:
        if self._hasParam:
            L = self.ParamGrid()
            x,y,z = L
            #print x,y,z
            #s = mlab.mesh(xb, yb, zb, scalars=real((Eplot_bulb*E*uxx).reshape(xb.shape)))
            # TODO: mayavi bug, opacity<1
            s = mlab.mesh(x, y, z, scalars=z, opacity=1.0)

        a,b = self.getBB()
        ll = (b+a)/2 - (extra_bb)*(b-a)/2
        rr = (b+a)/2 + (extra_bb)*(b-a)/2
        for i in range(0,200):
            x = np.random.uniform( ll, rr )
            cp,dist,bdy,other = self.cp(x)
            drawplot = False
            if bdy==1:
                if (np.random.random(1) < 1.0):
                    drawplot = True
                    col = 'g'
                    colt = (0.5,1,.2)
                    op = 0.9
            elif bdy==2:
                if (np.random.random(1) < 1.0):
                    drawplot = True
                    col = 'b'
                    colt = (.2,.5,1)
                    op = 0.9
            else:
                if ( (np.random.random(1) < 0.5) and (dist <= max((b-a)/10)) ):
                    drawplot = True
                    col = 'k'
                    colt = (0.5,0.5,0.5)
                    op = 0.3
            if drawplot:
                l = mlab.plot3d([x[0],cp[0]], [x[1],cp[1]], [x[2], cp[2]], color=colt, opacity=op)#, tube_radius=0.1)
        #mlab.title('3d viz test')
        mlab.show()

