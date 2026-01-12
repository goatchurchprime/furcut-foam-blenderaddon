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
    def __init__(self, layernumber, lengalong, cpt, rampdrop=0.0):
        self.layernumber = layernumber
        self.lengalong = lengalong
        self.cpt = cpt
        self.rampdrop = rampdrop
        self.arcengagements = [ ]
        
    def setcolumnarcengagements(self, clayers, rad, stockpt, ptprev):
        for i in range(self.layernumber, len(clayers)):
            layerz = clayers[i][0][0].z
            if (not self.arcengagements) or layerz > self.arcengagements[-1].z:
                self.arcengagements.append(ArcEngagement(P2.ConvertLZ(self.cpt), layerz, rad, stockpt))
                if ptprev:
                    self.arcengagements[-1].samearcslice(ptprev)
    

class CutPathColumns:
    def __init__(self, pcs, maxgap, mingap, rampgap, clayers, toolrad, stockpt):
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
                    lam = j*1.0/subdivs                    
                    mpt = Along(lam, prevcpt, pcs[i].cpt)
                    rampdrop = 0 if gap > rampgap else (pcs[i].cpt.z - mpt.z)
                    mpt = P3.ConvertCZ(mpt, pcs[i].cpt.z)
                    self.pcs.append(ContPointColumn(pcs[i].layernumber, -1, mpt, rampdrop))
            if gap > mingap:
                self.pcs.append(pcs[i])
        self.prevclearpc = None
        
        pcprev = None
        for pc in self.pcs:
            pc.setcolumnarcengagements(clayers, toolrad, stockpt, P2.ConvertLZ(pcprev.cpt) if pcprev else None)
            pcprev = pc

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

    def setxyrg(self):
        def mergerg(rg1, rg2):
            if rg1.hi == -10:
                return rg2
            if rg2.hi == -10:
                return rg1
            return I1(min(rg1.lo, rg2.lo), max(rg1.hi, rg2.hi))
        xrg, yrg = I1(0,-10), I1(0,-10)
        for j in range(len(self.pcs)):
            xrg = mergerg(xrg, self.pcs[j].arcengagements[0].xyrgpt(True))
            yrg = mergerg(yrg, self.pcs[j].arcengagements[0].xyrgpt(False))
        self.xrg, self.yrg = xrg, yrg

    def sliceoutarc(self, pt):
        bchange = False
        for pc in self.pcs:
            i = len(pc.arcengagements)-1
            while i >= 0 and pc.arcengagements[i].z >= pt.z:
                bchange = pc.arcengagements[i].samearcslice(P2.ConvertLZ(pt)) or bchange
                i -= 1
        return bchange

    def extractpath(self, withrampdrop):
        pth = [ ]
        nextprevclear = None
        nextprevclearj = 0
        for j in range(len(self.ipc)):
            pc = self.pcs[j]
            i = self.ipc[j]
            pthz = pc.arcengagements[i].z
            if withrampdrop:
                pthz -= pc.rampdrop
                if len(pth) != 0 and pthz < pth[-1].z:
                    pthz = pth[-1].z
            pth.append(P3.ConvertGZ(pc.arcengagements[i].pt, pthz))
            if i == 0 and (nextprevclear == None or j == nextprevclearj+1):
                nextprevclear = pc.arcengagements[0].pt
                nextprevclearj = j
            del pc.arcengagements[i:]
        if self.prevclearpc and len(pth) != 0:
            pth.insert(0, P3.ConvertGZ(self.prevclearpc, pth[0].z))
        self.prevclearpc = nextprevclear
        self.ipc.clear()
        for j in range(len(self.pcs)-1, -1, -1):
            if len(self.pcs[j].arcengagements) == 0:
                del self.pcs[j]
        return pth


class CutPathColumns2:
    def __init__(self, pcs, maxgap, mingap, rampgap, clayers, toolrad, stockpt):
        i0 = 0
        while not (pcs[i0] and pcs[i0].cpt):
            i0 += 1
        shaftrad = max(4*0.001, toolrad*0.667)
        print("ShaftRad", shaftrad)
        self.toolrad = toolrad
        self.pscs = [ ArcEngagement(P2.ConvertLZ(pcs[i0].cpt), pcs[i0].cpt.z, shaftrad, stockpt) ]
        

        for i in range(i0+1, len(pcs)):
            prevcpt = P3.ConvertGZ(self.pscs[-1].pt, self.pscs[-1].z)
            gap = (pcs[i].cpt - prevcpt).LenLZ()
            if gap > maxgap:
                subdivs = int(gap/maxgap) + 2
                for j in range(1, subdivs):
                    lam = j*1.0/subdivs                    
                    mpt = Along(lam, prevcpt, pcs[i].cpt)
                    rampdrop = 0 if gap > rampgap else (pcs[i].cpt.z - mpt.z)
                    mpt = P3.ConvertCZ(mpt, pcs[i].cpt.z)
                    self.pscs.append(ArcEngagement(P2.ConvertLZ(mpt), mpt.z, shaftrad, stockpt))
                    self.pscs[-1].rampdrop2 = rampdrop
            if gap > mingap:
                self.pscs.append(ArcEngagement(P2.ConvertLZ(pcs[i].cpt), pcs[i].cpt.z, shaftrad, stockpt))
        self.clearedpscs = -1

    def cutpasseng(self, maxengage):
        self.engagedpscs = self.clearedpscs + 1
        while self.engagedpscs < len(self.pscs) and self.pscs[self.engagedpscs].engagement() == 0.0:
            self.engagedpscs += 1

    def iscore2(self, pt):
        segi = self.engagedpscs - 1 - self.clearedpscs
        if segi == 0:
            return -1
        linkdist = (self.pscs[max(0, self.clearedpscs)].pt - pt).Len()
        return segi + max(0, 20-linkdist*1000)*2

    def setxyrg(self):
        def mergerg(rg1, rg2):
            if rg1.hi == -10:
                return rg2
            if rg2.hi == -10:
                return rg1
            return I1(min(rg1.lo, rg2.lo), max(rg1.hi, rg2.hi))
        xrg, yrg = I1(0,-10), I1(0,-10)
        for j in range(len(self.pscs)):
            xrg = mergerg(xrg, self.pscs[j].xyrgpt(True))
            yrg = mergerg(yrg, self.pscs[j].xyrgpt(False))
        self.xrg, self.yrg = xrg, yrg

    def sliceoutarc(self, pt):
        bchange = False
        for j in range(len(self.pscs)):
            if self.pscs[j].z >= pt.z:
                bchange = self.pscs[j].diffarcslice(P2.ConvertLZ(pt), self.toolrad) or bchange
        return bchange

    def extractpath(self, D):
        pth = [ ]
        for j in range(max(0, self.clearedpscs), self.engagedpscs):
            pth.append(P3.ConvertGZ(self.pscs[j].pt, self.pscs[j].z - self.pscs[j].rampdrop2))
            self.pscs[j].arcsegs.clear()
        self.clearedpscs = self.engagedpscs - 1
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
            #i0 = -1 # disable
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

    for i in range(len(pcsslices)):
        pcs = pcsslices[i]
        while len(pcs) >= 2 and pcs[0].cpt == None and pcs[1].cpt == None:
            print("del lead ", i, len(pcs))
            del pcs[0]

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
        maxengage = 80  # degrees
        for pcsT in pcsslicesT:
            cpc = CutPathColumns2(pcsT, maxgap=1.0*0.001, mingap=0.01*0.001, rampgap=4*0.001, clayers=clayers, toolrad=toolrad, stockpt=stockpt)
            cpc.cutpasseng(maxengage)
            cpc.setxyrg()
            pcsslices.append(cpc)
            
        slicestodraw = [ ]
        withrampdrop = True
        prevpt = P2(0,0)
        for t in range(320):
            bestcut = max(range(len(pcsslices)), key=lambda X:pcsslices[X].iscore2(prevpt))
            #bestcut = max(range(47-10, 47+10), key=lambda X:pcsslices[X].iscore)
            print("bestcut", bestcut, pcsslices[bestcut].iscore2(prevpt))
            cpc = pcsslices[bestcut]
            pth = cpc.extractpath(withrampdrop)
            if len(pth) <= 1:
                break
            cpc.cutpasseng(maxengage)
            cpc.setxyrg()
            prevpt = P2.ConvertLZ(pth[-1])

            slicestodraw.clear()
            slicestodraw.append(bestcut)

            pxrg = I1.AbsorbList(pt.x for pt in pth).Inflate(toolrad)
            pyrg = I1.AbsorbList(pt.y for pt in pth).Inflate(toolrad)
            print("pprg", pxrg, pyrg)
            for i in range(len(pcsslices)):
                cpc = pcsslices[i]
                #print(i, cpc.xrg, cpc.yrg)
                if (i != bestcut) and (cpc.xrg.lo < pxrg.hi and pxrg.lo < cpc.xrg.hi) and (cpc.yrg.lo < pyrg.hi and pyrg.lo < cpc.yrg.hi):
                    bchange = False
                    for pt in pth:
                        bchange = pcsslices[i].sliceoutarc(pt) or bchange
                    if bchange:
                        slicestodraw.append(i)
                        cpc.cutpasseng(maxengage)
                        cpc.setxyrg()
            print("slicestodraw", slicestodraw)
            #if len(cpc.pcs) == 0:
            #    del pcsslices[bestcut]

            #for i in range(len(pcsslices)):
            #pth = pcsslices[i].extractpath()

            mesh = bpy.data.meshes.new("dong%d" % t)
            mvertices = pth
            medges = [ (a, a+1)  for a in range(len(pth)-1) ]
            mesh.from_pydata(mvertices, medges, [])
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mobj.display_type = 'WIRE'
            uptoolpath.objects.link(mobj)

        #slicestodraw = range(len(pcsslices))
        #for i in slicestodraw:
        #    cpc = pcsslices[i]
        #    for ae in cpc.pscs:
        #        ae.makeengementmesh(uptoolpath)

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
        self.rampdrop2 = 0
    def halfplaneslice(self, normout, d):
        a = (d - P2.Dot(normout, self.pt))/self.rad
        if a <= -1:
            subarc(self.arcsegs, 0, 360)
        elif a < 1:
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
    def diffarcslice(self, apt, crad):
        vec = apt - self.pt
        vlensq = vec.Lensq()
        if vlensq == 0:
            assert (crad >= self.rad)
            self.arcsegs.clear()
            return True
        vlen = math.sqrt(vlensq)
        a = 0.5*(1 - ((crad*crad - self.rad*self.rad)/vlensq))
        ac = a*vlen/self.rad
        if ac < 1:
            if ac > -1:
                aa = math.degrees(math.acos(ac))
                an = vec.Arg()
                if an < 0:
                    an += 360
                return subarc(self.arcsegs, an-aa, an+aa)
            else:
                return subarc(self.arcsegs, 0, 360)
        return False
        
    def engagement(self):
        return sum((self.arcsegs[i+1]-self.arcsegs[i]  for i in range(0, len(self.arcsegs), 2)), 0)
    def xyrgptgen(self, bx):
        for i in range(len(self.arcsegs)):
            if bx:
                yield self.pt.u + math.cos(math.radians(self.arcsegs[i]))*self.rad
                if (i%2) == 0:
                    if self.arcsegs[i] < 180 < self.arcsegs[i+1]:
                        yield self.pt.u - self.rad
            else:
                yield self.pt.v + math.sin(math.radians(self.arcsegs[i]))*self.rad
                if (i%2) == 0:
                    if self.arcsegs[i] < 90 < self.arcsegs[i+1]:
                        yield self.pt.v + self.rad
                    if self.arcsegs[i] < 270 < self.arcsegs[i+1]:
                        yield self.pt.v - self.rad
                    
    def xyrgpt(self, bx):
        return I1.AbsorbList(self.xyrgptgen(bx)) if self.arcsegs else I1(0,-10)

    def makeengementmesh(self, collection, k=0, d=10):
        mvertices = [ ]
        medges = [ ]
        for i in range(0, len(self.arcsegs), 2):
            dhi = self.arcsegs[i+1]
            dlo = self.arcsegs[i]
            subdivs = int((dhi - dlo)/d)+1
            for j in range(0, subdivs+1):
                dd = Along(j*1.0/subdivs, dlo, dhi)
                mvertices.append(P3(math.cos(math.radians(dd))*self.rad + self.pt.u, math.sin(math.radians(dd))*self.rad + self.pt.v, self.z))
                if j != 0:
                    medges.append((len(mvertices)-2, len(mvertices)-1))
        if mvertices:
            mesh = bpy.data.meshes.new("zarc%d" % k)
            mesh.from_pydata(mvertices, medges, [])
            mobj = bpy.data.objects.new(mesh.name, mesh)
            mobj.display_type = 'WIRE'
            collection.objects.link(mobj)
