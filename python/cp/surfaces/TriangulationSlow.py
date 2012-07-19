"""
Closest Point Representation from a triangulation by brute force search.

TODO: have better approaches coded elsewhere

TODO: code from oo/cpFromTri.py, should be in class, separate file or what?  Triangle_ops as a module seems like a good idea
"""
from Surface import Surface
import numpy
from numpy import array as a
from numpy.linalg import norm

def ProjectOnSegment(c1in, c2in, c3in, p1, p2, p3, q1, q2, q3):
    """copied from steve's C code"""
    #double *c1, *c2, *c3;  <-- return values
    #double p1, p2, p3;
    #double q1, q2, q3;

    #double lamb, lamb_star;
    #double cmp1, cmp2, cmp3;
    #double qmp1, qmp2, qmp3;

    cmp1 = c1in-p1
    cmp2 = c2in-p2
    cmp3 = c3in-p3
    qmp1 =  q1-p1
    qmp2 =  q2-p2
    qmp3 =  q3-p3

    lamb = (cmp1*qmp1+cmp2*qmp2+cmp3*qmp3)/(qmp1*qmp1+qmp2*qmp2+qmp3*qmp3)
    lamb_star = max(0.,min(lamb,1.))

    #printf("lamb_star p %g %g %g %g\n",lamb_star, p1,p2,p3);

    c1out = p1+lamb_star*qmp1
    c2out = p2+lamb_star*qmp2
    c3out = p3+lamb_star*qmp3

    return (c1out, c2out, c3out)



def FindClosestPointToOneTri(a1,a2,a3, aface, vertex):
    """ copied from steve's C code """
    #double a1, a2, a3;
    #double *c1, *c2, *c3;   <-- mutable, pointers in C code, we will return them instead
    #int face_index;

    # int index_p, index_q, index_r;
    # double a11, a12, a22, b1, b2;
    # double i11, i12, i22, factor;
    # double q1, q2, q3;
    # double r1, r2, r3;
    # double lamb, mu;
    # double dd;
    # int i;

    # /* obtain the indices to the three vertices */
    #index_p = face[face_index].v1;
    #index_q = face[face_index].v2;
    #index_r = face[face_index].v3;
    index_p = aface[0]
    index_q = aface[1]
    index_r = aface[2]

    #/* translate so the p is at the origin */
    a1 -= vertex[index_p][0];
    a2 -= vertex[index_p][1];
    a3 -= vertex[index_p][2];
    q1 = vertex[index_q][0]-vertex[index_p][0];
    q2 = vertex[index_q][1]-vertex[index_p][1];
    q3 = vertex[index_q][2]-vertex[index_p][2];
    r1 = vertex[index_r][0]-vertex[index_p][0];
    r2 = vertex[index_r][1]-vertex[index_p][1];
    r3 = vertex[index_r][2]-vertex[index_p][2];

    #/* evaluate the various matrix entries */
    a11 = q1*q1+q2*q2+q3*q3;
    a12 = q1*r1+q2*r2+q3*r3;
    a22 = r1*r1+r2*r2+r3*r3;
    b1  = a1*q1+a2*q2+a3*q3;
    b2  = a1*r1+a2*r2+a3*r3;

    #/* find the inverse matrix and solve for lamb and mu */
    factor = 1.0/(a11*a22-a12*a12);
    i11 =  a22*factor;
    i12 = -a12*factor;
    i22 =  a11*factor;
    lamb = i11*b1+i12*b2;
    mu     = i12*b1+i22*b2;
    c1 = lamb*q1+mu*r1;  # Cptr
    c2 = lamb*q2+mu*r2;  # Cptr
    c3 = lamb*q3+mu*r3;  # Cptr

    if ((lamb<0) and (mu<0) and (lamb+mu<=1)):
        c1 = c2 = c3 = 0.0
    elif ((lamb>=0) and (mu<0) and (lamb+mu<=1)):
        # some acceleration possible using cython
        #(c1,c2,c3) = tri.C_ProjectOnSegment(c1,c2,c3,0.,0.,0.,q1,q2,q3);
        (c1,c2,c3) = ProjectOnSegment(c1,c2,c3,0.,0.,0.,q1,q2,q3);
    elif ((lamb>=0) and (mu<0) and (lamb+mu>1)):
        c1 = q1;
        c2 = q2;
        c3 = q3;
    elif ((lamb>=0) and (mu>=0) and (lamb+mu>1)):
        #(c1,c2,c3) = tri.C_ProjectOnSegment(c1,c2,c3,q1,q2,q3,r1,r2,r3);
        (c1,c2,c3) = ProjectOnSegment(c1,c2,c3,q1,q2,q3,r1,r2,r3);
    elif ((lamb<0) and (mu>=0) and (lamb+mu>1)):
        c1 = r1;
        c2 = r2;
        c3 = r3;
    elif ((lamb<0) and (mu>=0) and (lamb+mu<=1)):
        #(c1,c2,c3) = tri.C_ProjectOnSegment(c1,c2,c3,r1,r2,r3,0.,0.,0.);
        (c1,c2,c3) = ProjectOnSegment(c1,c2,c3,r1,r2,r3,0.,0.,0.);
    elif ((lamb>=0) and (mu>=0) and (lamb+mu<=1)):
        # /* do nothing */
        True
    else:
        print("non-enumerated case.\n");
        print("lamb mu %g %g\n", lamb, mu);
        raise NameError('case should not occur')

    #/* Calculate sqr distance( */
    dd  = (a1-c1)*(a1-c1)+(a2-c2)*(a2-c2)+(a3-c3)*(a3-c3);

# /*  DEBUGGING
# printf("lamb mu %g %g\n",lamb,mu);
# printf("q %g %g %g\n",q1,q2,q3);
# printf("r %g %g %g\n",r1,r2,r3);
# 	printf("%g %g %g %g\n",dd, *c1, *c2, *c3);
# printf("c %g %g %g\n",*c1,*c2,*c3);
# printf("a11 a12 a22 %g %g %g\n",a11,a12,a22);
# printf("c %g %g %g\n",*c1,*c2,*c3);

# printf("p %g %g %g\n",vertex[index_p][0],vertex[index_p][1],vertex[index_p][2]);
# printf("a %g %g %g\n",a1,a2,a3);
# printf("%g\n",dd);
# for (i=0;i<1000000;i++)
# {
# 	double w1,w2;
# 	w1=((double)(rand()%640001))/640000.;
# 	w2=((double)(rand()%640001))/640000.;
#         *c1 = w1*q1+(1-w1)*w2*r1;
#         *c2 = w1*q2+(1-w1)*w2*r2;
#         *c3 = w1*q3+(1-w1)*w2*r3;
#         dd  = (a1-*c1)*(a1-*c1)+(a2-*c2)*(a2-*c2)+(a3-*c3)*(a3-*c3);
# 	printf("%g %g %g %g\n",dd, *c1, *c2, *c3);
# }
# */
    #/* Shift everything back */
    c1 += vertex[index_p][0];
    c2 += vertex[index_p][1];
    c3 += vertex[index_p][2];
    return (dd,c1,c2,c3);



def FindClosestPointToTriSet(a1,a2,a3, Faces, Vertices):
    from numpy import inf
    #dd_min = 1.0e42  # just a big number
    dd_min = inf
    for F in Faces:
        (dd,t1,t2,t3) = FindClosestPointToOneTri(a1,a2,a3, F, Vertices)
        if (dd<dd_min):
            #minF = F
            dd_min = dd
            c1 = t1
            c2 = t2
            c3 = t3
    #print "minF = " + str(minF)
    return (dd_min,c1,c2,c3)


def loadPly(fname, fptype=numpy.float64):
    """
    load a (non-standard) ply file
    """
    from numpy import zeros
    f =  open(fname, 'r')

    # how many Vertices, Faces
    t = f.readline()
    t2 = t.split()
    numV = int(t2[0])
    numF = int(t2[1])
    print "number of vertices, faces = " + str((numV,numF))

    Vertices = numpy.zeros((numV, 3), dtype=fptype)
    #Vertices[:,0] = x
    #Vertices[:,1] = y
    #Vertices[:,2] = z
    for i in range(0, numV):
        t = f.readline().split()
        # TODO: suffers from float96 bug
        x = fptype(t[0])
        y = fptype(t[1])
        z = fptype(t[2])
        Vertices[i,0] = x
        Vertices[i,1] = y
        Vertices[i,2] = z
    Faces = numpy.zeros((numF,3), dtype=int)
    for i in range(0, numF):
        t = f.readline().split()
        # in case there is a quad hiding in there
        #if t != 'tri 3\r\n':
        #    print "help!!!"
        #t = f.readline()
        #t2 = t.split()
        # t[0] is 3
        Faces[i,0] = int(t[1])
        Faces[i,1] = int(t[2])
        Faces[i,2] = int(t[3])
    f.close
    return (Vertices,Faces)


class TriangulationSlow(Surface):
    def __init__(self):
        # todo: use a ply library?
        self.V, self.F = loadPly('bunny35947_input.ply')

    def closestPointToCartesian(self, xx, verbose=0):
        """ brue force search over all triangles """
        x,y,z=xx

        (dd,cpx,cpy,cpz) = FindClosestPointToTriSet(x,y,z, self.F, self.V)
        cp = a([cpx,cpy,cpz])
        dist = norm(xx - cp, 2)
        bdy = 0  # TODO: I guess!
        #others = dict(whichFace=tmin)
        others = {}
        if (verbose >=0):
            print "found cp: x,cp,dist=", xx,cp,dist
        return cp, dist, bdy, others

    cp = closestPointToCartesian
