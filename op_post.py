import bpy
import bmesh
import mathutils

print("reloading oppost ")
print("reloading opupcut ")
from .barmesh.basicgeo import I1, Partition1, P3, P2, Along

from bpy.types import Operator
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
)


def extractcurve(zslice, xoffset_mirrored):
    averts = [P3(xoffset_mirrored - v.co[0] if xoffset_mirrored else v.co[0], v.co[1], v.co[2])*1000  for v in zslice.data.vertices]
    edges = [tuple(e.vertices)  for e in zslice.data.edges]
    dedges = dict(edges)
    iv = iv0 = edges[0][0]
    cont = [ averts[iv0] ]
    while (iv := dedges.get(iv)) and iv != iv0:
        cont.append(averts[iv])
    return cont


class PostProc(bpy.types.Operator):
    bl_idname = "object.furcut_post"
    bl_label = "Furcut post"
    bl_options = {'REGISTER', 'UNDO'}

    mirrored: BoolProperty(
        name="mirror",
        description="mirror",
        default=False
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
        uptoolpath = bpy.data.collections[bpy.data.collections.find("uptoolpath")]

        zclear = 25.4+2
        ztop = 25.4
        nslab = 1

        stockpt = bpy.data.collections[bpy.data.collections.find("cncwork")].objects[1].data.vertices[2].co
        sdim = "%dx%d" % (int(stockpt[0]*1000+0.75), int(stockpt[1]*1000+0.75))
        if self.mirrored:
            sdim = sdim + "_mirrored"
        fname = "/home/julian/repositories/furcut/ncfiles/slab%d_%s.nc" % (nslab, sdim)

        fout = open(fname, "w")
        print(fname)

        fout.write("""%
N10 G90 G94 G17 G21
G55 ( workpiece offset )
N35 M9
N40 T1 M6
N45 S5000 M3
N55 M8 
G1 Z30F4000

""")

        topslide = 4  # this can be guided by the top rim engagement
        xoffset_mirrored = stockpt[0] if self.mirrored else 0
        for pth in uptoolpath.objects:
            vxs = extractcurve(pth, xoffset_mirrored)
            fout.write("X%.2f Y%.2f\n" % (vxs[0][0], vxs[0][1]))
            for p in vxs:
                fout.write("X%.2f Y%.2f Z%.2f\n" % (p[0], p[1], p[2]))
            if topslide:
                fout.write("Z%.3f\n" % ztop)
                v = P2.ZNorm(P2(p[0] - vxs[0][0], p[1] - vxs[0][1]))*topslide
                fout.write("X%.2f Y%.2f\n" % (p[0] + v.u, p[1] + v.v))
            fout.write("Z%.2f\n" % (zclear))
    
        fout.write("""X280Z40
N2670 M9
N2680 G90
N2695 M30
%
""")
        fout.close()

        return {'FINISHED'}

