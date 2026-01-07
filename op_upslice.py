import math
import bpy, bmesh
import mathutils
from mathutils import Vector
from bl_ext.blender_org.stl_format_legacy import blender_utils

#from .barmesh.basicgeo import I1, Partition1, P3, P2, Along
#import .barmesh, barmesh.implicitareaballoffset, barmesh.implicitareacyloffset
#from barmesh.barmeshslicer import BarMeshSlicer
#import barmesh.mainfunctions

#import sys
#sys.path.append("/home/julian/repositories/furcut/furcut-foam-blenderaddon/")
#from barmesh.basicgeo import I1, Partition1, P3, P2, Along

print("reloading opupslice ")
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
    pass
    return mcoll

def extractcurve(zslice):
    averts = [P3(*v.co)  for v in zslice.data.vertices]
    edges = [tuple(e.vertices)  for e in zslice.data.edges]
    dedges = dict(edges)
    iv = iv0 = edges[0][0]
    cont = [ averts[iv0] ]
    while (iv := dedges.get(iv)) and iv != iv0:
        cont.append(averts[iv])
    return cont

def convertconseq(cont):
    constseq = [ ]
    scont = cont[1:] + ([ cont[0] ] if cont[-1] != cont[0] else [ ]) 
    for cp0, cp1 in zip(cont, scont):
        constseq.append((cp0, cp1, (cp0 - cp1).Len() + (constseq[-1][2] if constseq else 0.0)))
    return constseq

def clengthalong(l, constseq):
    while l > constseq[-1][2]:
        l -= constseq[-1][2]
    for i in range(len(constseq)):
        if l <= constseq[i][2]:
            prevac = constseq[i-1][2] if i > 0 else 0.0
            lam = (l - prevac)/(constseq[i][2] - prevac)
            return Along(lam, constseq[i][0], constseq[i][1])
    assert False

def pointclosest(p, constseq, tanvec=None):
    assert (p != None)
    tanvecperp = None
    if tanvec:
        tanvecperp = P2.APerp(tanvec[1] - tanvec[0])
        tanvecval = P2.Dot(tanvecperp, tanvec[0])
    l = -1
    dsq = -1
    pc = None
    prevcl = 0.0
    for p0, p1, cl in constseq:
        if tanvec:
            tp0 = P2.DotLZ(tanvecperp, p0) - tanvecval
            tp1 = P2.DotLZ(tanvecperp, p1) - tanvecval
            if tp0 < 0 and tp1 < 0:
                continue
            if (tp0 < 0) != (tp1 < 0):
                mu = (0-tp0)/(tp1-tp0)
                tpm = Along(mu, p0, p1)
                if tp0 < 0:
                    tp0 = tpm
                else:
                    tp1 = tpm
        v = p1 - p0
        vsq = v.Lensq()
        if vsq != 0:
            v0 = p - p0
            lam = P3.Dot(p - p0, v)/vsq
            if lam < 1.0:
                if lam < 0.0:
                    ldsq = v0.Lensq()
                    lpc = p0
                else:
                    lpc = p0 + v*lam
                    ldsq = (lpc - p).Lensq()
                if pc == None or ldsq < dsq:
                    dsq = ldsq
                    l = Along(max(lam, 0.0), prevcl, cl)
                    pc = lpc
        prevcl = cl
    return pc, l


class ContPointColumn:
    def __init__(self, layernumber, lengalong, cpt):
        self.layernumber = layernumber
        self.lengalong = lengalong
        self.cpt = cpt
    
def trackpcsup(pcsprev, i0, l0, clayers, pcsslices):
    print("trackpcsup", len(pcsprev), len(clayers))
    pcs = [ ]
    for i in range(0, i0):
        pcs.append(ContPointColumn(i, pcsprev[i].lengalong if pcsprev else 0, None))
    cp0 = clengthalong(l0, clayers[i0])
    pcs.append(ContPointColumn(i0, l0, cp0))
    for i in range(i0+1, len(clayers)):
        clayer = clayers[i]
        cp, l = pointclosest(pcs[-1].cpt, clayer)
        while pcsslices and l < pcsslices[-1][i].lengalong and abs(l + clayer[-1][2] - pcsslices[-1][i].lengalong) < abs(l - pcsslices[-1][i].lengalong):
            l += clayer[-1][2]
        pcs.append(ContPointColumn(i, l, cp))
    return pcs


def sliceup(clayers, di, contourstepover):
    l0 = 0
    pcsslices = [ ]
    mincontourstepover = contourstepover*0.3
    #mincontourstepover = 0 # disable
    maxcontourstepover = contourstepover*1.5
    l0max = clayers[di][-1][2]
    while l0 < l0max:
        pcs = trackpcsup([], di, l0, clayers, pcsslices)
        pcsnextstack = [ pcs ]
        while pcsnextstack:
            pcsprev = pcsslices[-1] if pcsslices else None
            pcs = pcsnextstack.pop()
            i0 = -1
            if pcsprev:
                for i in range(di+1, len(clayers)-1):
                    if pcs[i].cpt and pcs[i].lengalong - pcsprev[i].lengalong > maxcontourstepover:
                        i0 = i
                        break
            if i0 != -1:
                divs = int((pcs[i0].lengalong - pcsprev[i0].lengalong)/contourstepover + 0.5)
                pcsnextstack.append(pcs)
                if divs > 2:
                    print("divs", i0, divs, (pcsprev[i0].lengalong, pcs[i0].lengalong))
                for j in range(divs-1, 0, -1):
                    li = Along(j/divs, pcsprev[i0].lengalong, pcs[i0].lengalong)
                    pcsi = trackpcsup(pcsprev, i0, li, clayers, pcsslices)
                    pcsnextstack.append(pcsi)
                continue
            pcsslices.append(pcs)
        l0 += contourstepover


    assert (di == 1)  # project back out to the outer contour curve
    for i in range(len(pcsslices)):
        pcs = pcsslices[i]
        pcsprev = pcsslices[i-1 if i else len(pcsslices)-1]
        pcsnext = pcsslices[i+1 if i!=len(pcsslices)-1 else 0]
        tanvec = (P2.ConvertLZ(ptlastpt(pcsprev)), P2.ConvertLZ(ptlastpt(pcsnext)))
        cpsideclearlayer, lsideclearlayer = pointclosest(pcs[di].cpt, clayers[0], tanvec)
        pcs[0] = ContPointColumn(0, lsideclearlayer, cpsideclearlayer)
    return pcsslices

def ptlastpt(pcs):
    i = len(pcs)-1
    while not pcs[i].cpt:
        i -= 1
    return pcs[i].cpt 

toolrad = (38.5/2)*0.001

class Upslice(bpy.types.Operator):
    bl_idname = "object.furcut_upslice"
    bl_label = "Furcut upslice"
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
        print("Bongo execute")
        zslicecollection = fetchcreatecollection("zslices", empty=False)
        stockpt = bpy.data.collections[bpy.data.collections.find("cncwork")].objects[1].data.vertices[2]
        uptoolpath = fetchcreatecollection("uptoolpath", empty=True)

        clayers = [ ]
        for zslice in zslicecollection.objects:
            cont = extractcurve(zslice)
            clayers.append(convertconseq(cont))

        stockpt = bpy.data.collections[bpy.data.collections.find("cncwork")].objects[1].data.vertices[2].co
        sideclearlayer = clayers[0]

        contourstepover = 6*0.001
        forcevertstepdist = 12*0.001
        pcsslices = sliceup(clayers, 1, contourstepover)
        
        pths = [ ]
        for i in range(len(pcsslices)):
            pth = [ pc.cpt  for pc in pcsslices[i]  if pc.cpt ]

            mesh = bpy.data.meshes.new("dong%d" % i)
            mvertices = pth
            medges = [ (a, a+1)  for a in range(len(pth)-1) ]
            mesh.from_pydata(mvertices, medges, [])
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mobj.display_type = 'WIRE'
            uptoolpath.objects.link(mobj)

            pths.append(pth)  


        return {'FINISHED'}

#lineslice
#circleslice
#subtractrange

# then we have a stack of these per segment along the line, including the stepin at the top
# the stacks start at the top but don't go all the way down
# Stack gets more full of material as we go down, so we look for the lowest pass in each stack that is within the engagement angle
# When we select a sequence stack we clear it and everything above
# Till we have passed along each.
# Update the adjacent stacks as we go along
# Choose the meatiest stack for next where it goes in the most for the least cusp engagement
# Engagement should be in direction of the motion
# Undercuts to be cut separately
#We look for a pass through them that picks the lowest point in each stack

def subarc(arcsegs, alo, ahi):
    if alo < 0:
        subarc(arcsegs, 360+alo, 360)
        alo = 0
    if ahi > 360:
        subarc(arcsegs, 0, ahi-360)
        ahi = 360
    i = len(arcsegs) - 2
    while i >= 0 and arcsegs[i+1] > alo:
        if arcsegs[i] < ahi:
            if alo > arcsegs[i]:
                if ahi < arcsegs[i+1]:
                    arcsegs.insert(i, arcsegs[i+1])
                    arcsegs.insert(i, arcsegs[i+1])
                    arcsegs[i+2] = ahi
                arcsegs[i+1] = alo
            elif ahi < arcsegs[i+1]:
                arcsegs[i] = ahi
            else:
                del arcsegs[i:i+2]
        i -= 2

class arcengagement:
    def __init__(self, pt, rad, stockpt):
        self.pt = pt
        self.rad = rad
        self.arcsegs = [ 0.0, 360.0 ]
        self.halfplaneslice(P2(-1,0), 0)
        self.halfplaneslice(P2(0,1), stockpt[1])
        self.halfplaneslice(P2(1,0), stockpt[0])
        self.halfplaneslice(P2(0,-1), 0)
    def halfplaneslice(self, normout, d):
        a = (d - P2.Dot(normout, self.pt))/self.rad
        if -1<a<1:
            aa = math.degrees(math.acos(a))
            an = normout.Arg()
            subarc(self.arcsegs, an-aa, an+aa)


