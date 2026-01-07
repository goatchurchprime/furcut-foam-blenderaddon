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
        self.arcengagements = [ ]
        
    def setcolumnarcengagements(self, clayers, rad, stockpt):
        for i in range(self.layernumber, len(clayers)):
            layerz = clayers[i][0][0].z
            if (not self.arcengagements) or layerz > self.arcengagements[-1].z:
                self.arcengagements.append(ArcEngagement(P2.ConvertLZ(self.cpt), layerz, rad, stockpt))
    

class CutPathColumns:
    def __init__(self, pcs, maxgap, mingap, clayers, toolrad, stockpt):
        i0 = 0
        while not (pcs[i0] and pcs[i0].cpt):
            i0 += 1
        self.pcs = [ pcs[i0] ]
        for i in range(i0+1, len(pcs)):
            prevcpt = self.pcs[-1].cpt
            gap = (pcs[i].cpt - prevcpt).LenLZ()
            if gap > maxgap:
                subdivs = int(gap/maxgap) + 2
                for j in range(1, subdivs):
                    self.pcs.append(ContPointColumn(pcs[i].layernumber, -1, P3.ConvertCZ(Along(j*1.0/subdivs, prevcpt, pcs[i].cpt), pcs[i].cpt.z)))
            if gap > mingap:
                self.pcs.append(pcs[i])
        self.prevclearpc = None
        
        for pc in self.pcs:
            pc.setcolumnarcengagements(clayers, toolrad, stockpt)
    
    def cutpasseng(self, maxengage):
        self.ipc = [ ]
        self.iscore = 0
        for j in range(len(self.pcs)):
            pc = self.pcs[j]
            i = len(pc.arcengagements)-1
            while i > 0 and (j == 0 or pc.arcengagements[i-1].z >= zcurrcut) and pc.arcengagements[i].engagement() < maxengage:
                i -= 1
            zcurrcut = pc.arcengagements[i].z
            self.iscore += len(pc.arcengagements) - i
            self.ipc.append(i)

    def sliceoutarc(self, pt):
        bchange = False
        for pc in self.pcs:
            i = len(pc.arcengagements)-1
            while i >= 0 and pc.arcengagements[i].z >= pt.z:
                bchange = pc.arcengagements[i].samearcslice(P2.ConvertLZ(pt)) or bchange
        return bchange

    def extractpath(self):
        pth = [ ]
        nextprevclear = None
        for j in range(len(self.ipc)):
            pc = self.pcs[j]
            i = self.ipc[j]
            pth.append(P3.ConvertGZ(pc.arcengagements[i].pt, pc.arcengagements[i].z))
            nextprevclear = pc.arcengagements[0].pt if i == 0 else None
            del pc.arcengagements[i:]
        if self.prevclearpc:
            pth.insert(0, P3.ConvertGZ(self.prevclearpc, pth[0].z))
        self.prevclearpc = nextprevclear
        self.ipc.clear()
        for j in range(len(self.pcs)-1, -1, -1):
            if len(self.pcs[j].arcengagements) == 0:
                del self.pcs[j]
        return pth
    
def trackpcsup(pcsprev, i0, l0, clayers, pcsslices):
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
            i0 = -1 # disable
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
        if pcs[di].cpt != None:
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
        pcsslicesT = sliceup(clayers, 1, contourstepover)
        pcsslices = [ ]
        maxengage = 70  # degrees
        for pcsT in pcsslicesT:
            cpc = CutPathColumns(pcsT, maxgap=1.0*0.001, mingap=0.01*0.001, clayers=clayers, toolrad=toolrad, stockpt=stockpt)
            cpc.cutpasseng(maxengage)
            pcsslices.append(cpc)
            
        for t in range(3):
            bestcut = max(range(len(pcsslices)), key=lambda X:pcsslices[X].iscore)
            print("bestcut", bestcut, pcsslices[bestcut].iscore)
            cpc = pcsslices[bestcut]
            #print("sss ", [ len(pc.arcengagements)  for pc in cpc.pcs ])
            pth = cpc.extractpath()
            #print("sssp ", [ len(pc.arcengagements)  for pc in cpc.pcs ])
            cpc.cutpasseng(maxengage)

            for i in range(len(pcsslices)):
                cpc = pcsslices[i]
                if i != bestcut:
                    bchange = False
                    for pt in pth:
                        bchange = pcsslices[i].sliceoutarc(pt) or bchange
                    if bchange:
                        cpc.cutpasseng(maxengage)

            #for i in range(len(pcsslices)):
            #pth = pcsslices[i].extractpath()

            mesh = bpy.data.meshes.new("dong%d" % t)
            mvertices = pth
            medges = [ (a, a+1)  for a in range(len(pth)-1) ]
            mesh.from_pydata(mvertices, medges, [])
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mobj.display_type = 'WIRE'
            uptoolpath.objects.link(mobj)

        return {'FINISHED'}

# Choose the meatiest stack for next where it goes in the most for the least cusp engagement
# Engagement should be in direction of the motion
# Undercuts to be cut separately
#We look for a pass through them that picks the lowest point in each stack

def subarc(arcsegs, alo, ahi):
    #print("subarc", alo, ahi)
    res = False
    if alo < 0:
        res = subarc(arcsegs, 360+alo, 360) or res
        alo = 0
    if ahi > 360:
        res = subarc(arcsegs, 0, ahi-360) or res
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
                res = True
            elif ahi < arcsegs[i+1]:
                arcsegs[i] = ahi
                res = True
            else:
                del arcsegs[i:i+2]
                res = True
        i -= 2
    return res

class ArcEngagement:
    def __init__(self, pt, z, rad, stockpt):
        self.pt = pt
        self.z = z
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
            if an < 0:
                an += 360
            subarc(self.arcsegs, an-aa, an+aa)
    def samearcslice(self, apt):
        vec = apt - self.pt
        a = vec.Len()*0.5/self.rad
        if a < 1:
            aa = math.degrees(math.acos(a))
            an = vec.Arg()
            if an < 0:
                an += 360
            return subarc(self.arcsegs, an-aa, an+aa)
        return False
    def engagement(self):
        return sum((self.arcsegs[i+1]-self.arcsegs[i]  for i in range(0, len(self.arcsegs), 2)), 0)

