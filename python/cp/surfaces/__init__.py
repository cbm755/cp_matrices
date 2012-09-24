"""
TODO: can we load all the various classes automatically here?
"""
from warnings import warn
# Surface is base classes only, don't import directly
from Hemisphere import Hemisphere
from Sphere import Sphere
from Circle import Circle
from CPBar import CPBar
from Mesh import Mesh
try:
    from triangulation_fast import FindClosestPointToTriSet
except ImportError:
    warn('fast triangulation code not compiled')

# TODO:
# treatment of __all__ for "from closestpoint import *"
