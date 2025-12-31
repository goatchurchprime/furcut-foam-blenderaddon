from .trianglebarmesh import TriangleBarMesh, TriangleBar
from .triangleboxing import MakeTriangleBoxing
from .trianglezcut import TriZCut

# doesn't work for reloading this
#import importlib, sys
#import barmesh.tribarmes.trianglebarmesh
#importlib.reload(barmesh.tribarmes.trianglebarmesh)
#TriangleBarMesh = barmesh.tribarmes.trianglebarmesh.TriangleBarMesh

bnumpyexists = True
try:  import numpy
except ImportError:  bnumpyexists = False

if bnumpyexists:
    from .ntrianglebarmesh import NTriangleBarMesh
else:
    NTriangleBarMesh = None
    
    
