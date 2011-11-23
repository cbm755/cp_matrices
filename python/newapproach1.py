import numpy as np


x1d,dx = np.linspace(-2,2,11,retstep=True)
y1d = x1d.copy()   # otherwise, its a pointer

[xx,yy] = np.meshgrid(x1d,y1d);

import closestpoint
from closestpoint import coordinate_transform as tf
q = closestpoint.Hemisphere()

# wrap the object in CPBar, for accurately imposing boundary
# conditions on shapes with boundaries (like a hemisphere)
#q = closestpoint.CPBar(parent=q)

cpfun = q.cp

th, r = tf.cart2pol(xx, yy)
cpx, cpy = tf.pol2cart(th, 1.1)
#cpx, cpy = tf.pol2cart(th, self._radius)


