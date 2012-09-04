"""
TODO: can we load all the various classes automatically here?
"""
# Surface is base classes only, don't import directly
from Hemisphere import Hemisphere
from Sphere import Sphere
from Circle import Circle
from Multi import Multi
from CPBar import CPBar
from Manipulate import Translate  # Doesn't seem very useful
from ParamCurve import ParamCurve
from SplineCurve import SplineCurve
from CircularArc import CircularArc
from Cylinder import Cylinder
from triangulation_fast import FindClosestPointToTriSet

# TODO:
# treatment of __all__ for "from closestpoint import *"
