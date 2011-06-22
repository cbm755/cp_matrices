"""
TODO: can we load all the various classes automatically here?
"""
#import Hemisphere
# ClosestPoint is base classes only, don't import directly
from Hemisphere import *
from Sphere import *
from Circle import *
from Multi import *
from CPBar import *
from Manipulate import *
from ParamCurve import *
from SplineCurve import *
from CircularArc import *
from Cylinder import *
from TriangulationSlow import *

#reload Sphere


# TODO:
# treatment of __all__ for "from closestpoint import *"
