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

def placejobmeshobject(collection, name, mesh):
    print("in placejobmeshobject")
    i = collection.objects.find(name)
    if i == -1:
        mobj = bpy.data.objects.new(name, mesh)
        collection.objects.link(mobj)
    else:
        mobj = collection.objects.objects[i]
        mobj.data = mesh

    mobj.rotation_mode = 'XYZ'
    facemax = max((p.area, p)  for p in mobj.data.polygons)[1]
    if (facemax.normal - Vector((-1,0,0))).length < 1e-6:
        mobj.rotation_euler = (0,math.radians(-90),0)
    elif (facemax.normal - Vector((1,0,0))).length < 1e-6:
        mobj.rotation_euler = (0,math.radians(90),0)
    else:
        print("Unknown normal direction ", facemax.normal)

    print("matrix before loc", mobj.matrix_basis)
    mobj.location = -facemax.center
    mobj.location.rotate(mobj.rotation_euler)
    #mobj.rotation_euler.rotate(mathutils.Euler((0,0,math.radians(90))))
    print("matrix after loc", mobj.matrix_basis)
    print("location", mobj.location, mobj.rotation_euler)

    # apply the transform (could be used to find the range and minimize it)
    matrix = mobj.matrix_basis.copy()
    print("apply matrix", matrix)
    for vert in mobj.data.vertices:
        vert.co = matrix @ vert.co

    mobj.matrix_basis.identity()
    return mobj


class Upcut(bpy.types.Operator):
    bl_idname = "object.furcut_upcut"
    bl_label = "Furcut test"
    bl_options = {'REGISTER', 'UNDO'}

    remove_poles_beforehand: BoolProperty(
        name="yipyip",
        description="yipee",
        default=True
    )

    min_length: FloatProperty(
        name="Min Edge Length",
        description="Edges below this length will disappear",
        default=0.02,
        min=0,
        soft_max=0.5
    )

    def execute(self, context):
        print("Bingo execute")
        obj = bpy.context.selected_objects[0]
        mesh = bpy.data.meshes.new_from_object(obj.evaluated_get(bpy.context.evaluated_depsgraph_get()))
        cncworkcollection = fetchcreatecollection("cncwork", empty=True)
        print("calling placejobmeshobject")
        mobj = placejobmeshobject(cncworkcollection, "jobmesh", mesh)

        print("mobj verts", len(mobj.data.vertices))
        tris = blender_utils.faces_from_mesh(mobj, mathutils.Matrix.Diagonal((1000,1000,1000,0)))
        tbarmesh = TriangleBarMesh()
        tbarmesh.BuildTriangleBarmesh(tris)
        print("nodes", len(tbarmesh.nodes))
        print("xrg", tbarmesh.xlo, tbarmesh.xhi)
        print("yrg", tbarmesh.ylo, tbarmesh.yhi)
        print("zrg", tbarmesh.zlo, tbarmesh.zhi)
        zss = sliceit(tbarmesh)
        zslicecollection = fetchcreatecollection("zslices")
        for i, zs in enumerate(zss):
            mesh = bpy.data.meshes.new("ding%d" % i)
            mvertices = zs
            medges = [ (a, a+1)  for a in range(len(zs)-1) ]
            mesh.from_pydata(mvertices, medges, [])
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mobj.display_type = 'WIRE'
            zslicecollection.objects.link(mobj)

        return {'FINISHED'}


discrad = 29.62/2
discheight = 1.6
zlo = 0
zhi = 25.4
cgoffsets = [ (1,0,0,0,False), (2.1,8,7.8,4,False), (3.2,16,16,8,True), 
              (4,23,23,8,True), (5.1,28,28,8,True)
               ]
zstep, hthickness, vthickness, contourext, trimboundary = cgoffsets[0]
zstep = 5

def sliceit(tbarmesh):
    zsteps = int((zhi - zlo)/zstep + 0.5)
    zlevels = [ Along(zi/zsteps, zlo, zhi)  for zi in range(zsteps+1) ]
    tdiscrad = discrad + hthickness
    contouroffset = discrad + contourext

    iaoffset = barmesh.implicitareacyloffset.ImplicitAreaCylOffset(tbarmesh)
    #iaoffset = barmesh.implicitareaballoffset.ImplicitAreaBallOffset(tbarmesh)
    rex = tdiscrad + 2.5
    xpart = Partition1(tbarmesh.xlo-rex, tbarmesh.xhi+rex, 49)
    ypart = Partition1(tbarmesh.ylo-rex, tbarmesh.yhi+rex, 37)

    res = [ ]
    for z in zlevels:
        print(z)
        iaoffset.SetCylZrg(max(zlo, z - vthickness), min(zhi, z + discheight + vthickness))

        bm = barmesh.barmesh.BarMesh()
        bm.BuildRectBarMesh(xpart, ypart, z)
        rd2 = max(xpart.vs[1]-xpart.vs[0], ypart.vs[1]-ypart.vs[0], tdiscrad*1.5) + 0.1
        bms = barmesh.barmeshslicer.BarMeshSlicer(bm, iaoffset, rd=tdiscrad, rd2=rd2, contourdotdiff=0.95, contourdelta=0.05, lamendgap=0.001)

        #bms.initializecutsanddistances() # quick version
        bms.fullmakeslice()

        #contsN, topbars = barmesh.mainfunctions.BarMeshContoursN(bm, barmesh.barmesh.PZ_BEYOND_R)
        contours, contbars = barmesh.mainfunctions.BarMeshContoursF(bm, barmesh.barmesh.PZ_BEYOND_R)

        contnest = barmesh.mainfunctions.NestContours(contbars, barmesh.barmesh.PZ_BEYOND_R)
        mconts = dict((topbar.midcontournumber, contours)  for cont, topbar in zip(contours, contbars))
        outerconts = [mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if izone == barmesh.barmesh.PZ_BEYOND_R and outxn == -1]
        innerconts = [mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if not (izone == barmesh.barmesh.PZ_BEYOND_R and outxn == -1)]
        print(contnest)  # contnest not working yet, so trimming just first contour out of the list thing
        for conts in outerconts:
            for cont in conts[:1]:
                res.append([p*0.001  for p in cont])

    return res
