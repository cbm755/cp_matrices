from numpy.testing import decorators

##@decorators.knownfailureif(True, "Fails because closespoint doesn't exist.")
@decorators.skipif(True, "Fails because closespoint doesn't exist.")
def test_cp():
    """Don't really know what this is suppposed to test."""
    import closestpoint as cp

    #q = cp.Circle(a([10,12]), 6)
    #q2 = cp.Semicircle(a([10,12]), 6)
    #q2.viztest()

    #q3 = cp.Hemisphere(a([1,1,-3]), 6)
    #q3.viztest()


    #q4 = cp.Sphere(a([1,1,-3]), 6)
    #q4.viztest()


    q = cp.Semicircle()
    #q.viztest()
    q2 = cp.CPBar(q)
    q2.viztest()

    q = cp.Hemisphere()
    #q.viztest()
    q3 = cp.CPBar(q)
    q3.viztest()
