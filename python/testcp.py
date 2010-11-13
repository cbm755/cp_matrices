from numpy import array as a
import closestpoint as cp

q = cp.Circle(a([10,12]), 6)

q2 = cp.Semicircle(a([10,12]), 6)

#q2.viztest()


q3 = cp.Hemisphere(a([1,1,-3]), 6)

q3.viztest()


q4 = cp.Sphere(a([1,1,-3]), 6)

q4.viztest()
