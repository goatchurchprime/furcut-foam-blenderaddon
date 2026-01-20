import math
import bpy, bmesh
import mathutils
from mathutils import Vector
from bl_ext.blender_org.stl_format_legacy import blender_utils

#from .barmesh.basicgeo import I1, Partition1, P3, P2, Along
#import .barmesh, barmesh.implicitareaballoffset, barmesh.implicitareacyloffset
#from barmesh.barmeshslicer import BarMeshSlicer
#import barmesh.mainfunctions

print("reloading opupcut ")
from . import barmesh
from .barmesh.basicgeo import I1, Partition1, P3, P2, Along
from .barmesh.tribarmes import TriangleBarMesh, TriangleBar, MakeTriangleBoxing

from bpy.types import Operator
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
)

def fetchcreatecollection(name, empty=True):
    fi = bpy.data.collections.find(name)
    if fi == -1:
        mcoll = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(mcoll)
    else:
        mcoll = bpy.data.collections[fi]
        if empty:
            for o in mcoll.objects:
                mcoll.objects.unlink(o)
                bpy.data.objects.remove(o, do_unlink=True)
    return mcoll

def xyrange(points):
    return [ min(v[0] for v in points), max(v[0] for v in points),
             min(v[1] for v in points), max(v[1] for v in points) ]

def getmeshwithname(collection, name, mesh):
    i = collection.objects.find(name)
    if i == -1:
        mobj = bpy.data.objects.new(name, mesh)
        collection.objects.link(mobj)
    else:
        mobj = collection.objects.objects[i]
        mobj.data = mesh
    return mobj
    
def placejobmeshobject(collection, mesh, hclear, flipz):
    print("in placejobmeshobject")
    mobj = getmeshwithname(collection, "jobmesh", mesh)
    mobj.location = (0,0,0)
    mobj.rotation_mode = 'XYZ'
    facemax = max(mobj.data.polygons, key=lambda X: X.area)
    print("Normal direction ", facemax.normal, facemax.center)
    if (facemax.normal - Vector((-1,0,0))).length < 1e-6:
        mobj.rotation_euler = (0,math.radians(-90),0)
    elif (facemax.normal - Vector((1,0,0))).length < 1e-6:
        mobj.rotation_euler = (0,math.radians(90),0)
    elif (facemax.normal - Vector((0,0,-1))).length < 1e-6:
        mobj.rotation_euler = (0,0,0)
    elif (facemax.normal - Vector((0,0,1))).length < 1e-6:
        mobj.rotation_euler = (0,math.radians(180),0)
    elif (facemax.normal - Vector((0,1,0))).length < 1e-6:
        mobj.rotation_euler = (math.radians(-90),0,0)
    elif (facemax.normal - Vector((0,-1,0))).length < 1e-6:
        mobj.rotation_euler = (math.radians(90),0,0)
    else:
        print("*** Unknown normal direction ", facemax.normal)

    if flipz:
        mobj.rotation_euler.y += math.radians(180)

    xyareas = [ ]
    for d in range(0,90,5):
        mobj.rotation_euler.z = math.radians(d)
        matrix = mobj.matrix_basis.copy()
        xyareas.append((d, xyrange([ (matrix @ vert.co)  for vert in mobj.data.vertices ])))
    X = min(xyareas, key=lambda X: ((X[1][1]-X[1][0])*(X[1][3]-X[1][2])))
    rot90 = ((X[1][1]-X[1][0]) < (X[1][3]-X[1][2]))
    mobj.rotation_euler.z = math.radians(X[0] + (90 if rot90 else 0))
    print("matrix before loc", rot90, mobj.rotation_euler)
    print("X", X)
    mobj.location = -facemax.center
    mobj.location.rotate(mobj.rotation_euler)
    matrix = mobj.matrix_basis.copy()
    lrg = xyrange([ (matrix @ vert.co)  for vert in mobj.data.vertices ])
    print(mobj.location, (-X[1][2],-X[1][3],0) if rot90 else (-X[1][0],X[1][2],0))
    mobj.location.x += hclear - lrg[0]
    mobj.location.y += hclear - lrg[2]

    if flipz:
        mobj.location.z += 25.4*0.001 # hardcoded thickness

    #mobj.location.rotate(mobj.rotation_euler)
    #mobj.rotation_euler.rotate(mathutils.Euler((0,0,math.radians(90))))
    print("location", mobj.location, mobj.rotation_euler)

    # apply the transform (could be used to find the range and minimize it)
    matrix = mobj.matrix_basis.copy()
    print("apply matrix", matrix)
    for vert in mobj.data.vertices:
        vert.co = matrix @ vert.co
    mobj.matrix_basis.identity()

    stockmesh = bpy.data.meshes.new("stock")
    sobj = getmeshwithname(collection, "stockmesh", stockmesh)
    scorn = (hclear*2 + lrg[1]-lrg[0], hclear*2 + lrg[3]-lrg[2])
    print("Stock dimensions width", scorn[0]*1000, "height", scorn[1]*1000)
    svertices = [ (0,0,0), (scorn[0],0,0), (scorn[0],scorn[1],0), (0,scorn[1],0) ]
    sedges = [ (0,1), (1,2), (2,3), (3,0) ]
    stockmesh.from_pydata(svertices, sedges, [])
    sobj.matrix_basis.identity()
    sobj.display_type = 'WIRE'

class UpPlace(bpy.types.Operator):
    bl_idname = "object.furcut_placer"
    bl_label = "Furcut placer"
    bl_options = {'REGISTER', 'UNDO'}

    horiz_clearance: FloatProperty(
        name="HClearance",
        description="Stock allowance around the part",
        default=0.005,
        min=0,
        soft_max=0.1
    )
    
    flipz: BoolProperty(
        name="flipz",
        description="flip over",
        default=False
    )

    #def invoke(self, context, event):
    #    wm = context.window_manager
    #    return wm.invoke_props_dialog(self, width=250)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Places selected part on origin.", icon='INFO')
        layout.separator()
        layout.row()
        layout.row()
        row = layout.row()
        row.prop(self, "horiz_clearance")
        row = layout.row()
        row.prop(self, "flipz")

    def execute(self, context):
        print("Bingo execute")
        obj = bpy.context.selected_objects[0]
        mesh = bpy.data.meshes.new_from_object(obj.evaluated_get(bpy.context.evaluated_depsgraph_get()))
        cncworkcollection = fetchcreatecollection("cncwork", empty=True)
        print("calling placejobmeshobject")
        mobj = placejobmeshobject(cncworkcollection, mesh, self.horiz_clearance, self.flipz)
        return {'FINISHED'}
