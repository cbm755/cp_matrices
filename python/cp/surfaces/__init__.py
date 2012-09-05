"""
TODO: can we load all the various classes automatically here?
"""
# Surface is base classes only, don't import directly
from Hemisphere import Hemisphere
from Sphere import Sphere
from Circle import Circle
from CPBar import CPBar
from Mesh import Mesh
from triangulation_fast import FindClosestPointToTriSet

# TODO:
# treatment of __all__ for "from closestpoint import *"
