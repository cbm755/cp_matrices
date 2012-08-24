import numpy as np
cimport numpy as np
cimport cython


@cython.cdivision(True)
cdef inline void ProjectOnSegment(double * c1, double * c2, double * c3, double p1, double p2, double p3, double q1, double q2, double q3):
    """copied from steve's C code"""
    #double *c1, *c2, *c3  <-- return values
    #double p1, p2, p3
    #double q1, q2, q3

    #double lamb, lamb_star
    #double cmp1, cmp2, cmp3
    #double qmp1, qmp2, qmp3

    cdef double cmp1 = c1[0]-p1
    cdef double cmp2 = c2[0]-p2
    cdef double cmp3 = c3[0]-p3
    cdef double qmp1 =  q1-p1
    cdef double qmp2 =  q2-p2
    cdef double qmp3 =  q3-p3

    cdef double lamb = (cmp1*qmp1+cmp2*qmp2+cmp3*qmp3)/(qmp1*qmp1+qmp2*qmp2+qmp3*qmp3)
    cdef double lamb_star = max(0.,min(lamb,1.))

    #printf("lamb_star p %g %g %g %g\n",lamb_star, p1,p2,p3)

    c1[0] = p1+lamb_star*qmp1
    c2[0] = p2+lamb_star*qmp2
    c3[0] = p3+lamb_star*qmp3

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void FindClosestPointToOneTri(double a1, double a2, double a3, int[::1] aface, double[:,::1] vertex, double * dd_, double * t1, double * t2, double * t3):
    """ copied from steve's C code """
    #double a1, a2, a3
    #double *c1, *c2, *c3   <-- mutable, pointers in C code, we will return them instead
    #int face_index

    # int index_p, index_q, index_r
    # double a11, a12, a22, b1, b2
    # double i11, i12, i22, factor
    # double q1, q2, q3
    # double r1, r2, r3
    # double lamb, mu
    # double dd
    # int i

    # /* obtain the indices to the three vertices */
    #index_p = face[face_index].v1
    #index_q = face[face_index].v2
    #index_r = face[face_index].v3
    cdef int index_p = aface[0]
    cdef int index_q = aface[1]
    cdef int index_r = aface[2]

    #/* translate so the p is at the origin */
    a1 -= vertex[index_p, 0]
    a2 -= vertex[index_p, 1]
    a3 -= vertex[index_p, 2]
    cdef double q1 = vertex[index_q, 0]-vertex[index_p, 0]
    cdef double q2 = vertex[index_q, 1]-vertex[index_p, 1]
    cdef double q3 = vertex[index_q, 2]-vertex[index_p, 2]
    cdef double r1 = vertex[index_r, 0]-vertex[index_p, 0]
    cdef double r2 = vertex[index_r, 1]-vertex[index_p, 1]
    cdef double r3 = vertex[index_r, 2]-vertex[index_p, 2]

    #/* evaluate the various matrix entries */
    cdef double a11 = q1*q1+q2*q2+q3*q3
    cdef double a12 = q1*r1+q2*r2+q3*r3
    cdef double a22 = r1*r1+r2*r2+r3*r3
    cdef double b1  = a1*q1+a2*q2+a3*q3
    cdef double b2  = a1*r1+a2*r2+a3*r3

    #/* find the inverse matrix and solve for lamb and mu */
    cdef double factor = 1.0/(a11*a22-a12*a12)
    cdef double i11 =  a22*factor
    cdef double i12 = -a12*factor
    cdef double i22 =  a11*factor
    cdef double lamb = i11*b1+i12*b2
    cdef double mu = i12*b1+i22*b2
    cdef double c1 = lamb*q1+mu*r1  # Cptr
    cdef double c2 = lamb*q2+mu*r2  # Cptr
    cdef double c3 = lamb*q3+mu*r3  # Cptr

    if ((lamb<0) and (mu<0) and (lamb+mu<=1)):
        c1 = c2 = c3 = 0.0
    elif ((lamb>=0) and (mu<0) and (lamb+mu<=1)):
        ProjectOnSegment(&c1,&c2,&c3,0.,0.,0.,q1,q2,q3)
    elif ((lamb>=0) and (mu<0) and (lamb+mu>1)):
        c1 = q1
        c2 = q2
        c3 = q3
    elif ((lamb>=0) and (mu>=0) and (lamb+mu>1)):
        ProjectOnSegment(&c1,&c2,&c3,q1,q2,q3,r1,r2,r3)
    elif ((lamb<0) and (mu>=0) and (lamb+mu>1)):
        c1 = r1
        c2 = r2
        c3 = r3
    elif ((lamb<0) and (mu>=0) and (lamb+mu<=1)):
        ProjectOnSegment(&c1,&c2,&c3,r1,r2,r3,0.,0.,0.)
    elif ((lamb>=0) and (mu>=0) and (lamb+mu<=1)):
        # /* do nothing */
        pass
#    else:
#        print("non-enumerated case.\n")
#        print("lamb mu %g %g\n", lamb, mu)
#        raise NameError('case should not occur')

    #/* Calculate sqr distance( */
    cdef double dd  = (a1-c1)*(a1-c1)+(a2-c2)*(a2-c2)+(a3-c3)*(a3-c3)

# /*  DEBUGGING
# printf("lamb mu %g %g\n",lamb,mu)
# printf("q %g %g %g\n",q1,q2,q3)
# printf("r %g %g %g\n",r1,r2,r3)
# 	printf("%g %g %g %g\n",dd, *c1, *c2, *c3)
# printf("c %g %g %g\n",*c1,*c2,*c3)
# printf("a11 a12 a22 %g %g %g\n",a11,a12,a22)
# printf("c %g %g %g\n",*c1,*c2,*c3)

# printf("p %g %g %g\n",vertex[index_p][0],vertex[index_p][1],vertex[index_p][2])
# printf("a %g %g %g\n",a1,a2,a3)
# printf("%g\n",dd)
# for (i=0i<1000000i++)
# {
# 	double w1,w2
# 	w1=((double)(rand()%640001))/640000.
# 	w2=((double)(rand()%640001))/640000.
#         *c1 = w1*q1+(1-w1)*w2*r1
#         *c2 = w1*q2+(1-w1)*w2*r2
#         *c3 = w1*q3+(1-w1)*w2*r3
#         dd  = (a1-*c1)*(a1-*c1)+(a2-*c2)*(a2-*c2)+(a3-*c3)*(a3-*c3)
# 	printf("%g %g %g %g\n",dd, *c1, *c2, *c3)
# }
# */
    #/* Shift everything back */
    c1 += vertex[index_p, 0]
    c2 += vertex[index_p, 1]
    c3 += vertex[index_p, 2]
    dd_[0] = dd
    t1[0] = c1
    t2[0] = c2
    t3[0] = c3

def FindClosestPointToTriSet(double a1, double a2, double a3, int[:, ::1] Faces, double[:, ::1] Vertices):
    cdef double dd_min = np.inf
    cdef double dd, c1, c2, c3, t1, t2, t3
    cdef Py_ssize_t i, n = len(Faces)
    cdef int[::1] F
    for i in range(n):
        F = Faces[i, :]
        FindClosestPointToOneTri(a1, a2, a3, F, Vertices, &dd, &t1, &t2, &t3)
        if dd < dd_min:
            dd_min = dd
            c1 = t1
            c2 = t2
            c3 = t3
    return dd_min, c1, c2, c3
