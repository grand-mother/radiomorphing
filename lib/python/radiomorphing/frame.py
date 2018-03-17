# Local frame transforms for pulse shape computations.
import numpy

def get_rotation(zen, az, phigeo, bfieldangle):
    """Utility function for getting the rotation matrix between frames
    """
    s = numpy.sin(bfieldangle)
    B = numpy.array([numpy.cos(phigeo) * s, numpy.sin(phigeo) * s,
                     numpy.cos(bfieldangle)])
    
    #print "bfield vector ", B
    
    
    s = numpy.sin(zen)
    #v = numpy.array([-numpy.cos(az) * s, -numpy.sin(az) * s, -numpy.cos(zen)])
    v = numpy.array([numpy.cos(az) * s, numpy.sin(az) * s, numpy.cos(zen)])


    vxB = numpy.cross(v, B)
    vxB /= numpy.linalg.norm(vxB)
    vxvxB = numpy.cross(v, vxB)
    vxvxB /= numpy.linalg.norm(vxvxB)
    
    #print "v ", v, " vxB ", vxB , " vxvxB ",vxvxB

       
    return numpy.array((v, vxB, vxvxB))

def UVWGetter(cx, cy, cz, zen, az, phigeo, bfieldangle):
    """Closure for getting coordinates in the shower frame.
    """
    R = get_rotation(zen, az, phigeo, bfieldangle)
    origin = numpy.array((cx, cy, cz))

    def GetUVW(pos):
       return numpy.dot(R, pos - origin)
    return GetUVW

def XYZGetter(cx, cy, cz, zen, az, phigeo, bfieldangle):
    """Closure for getting back to the main frame
    """
    Rt = get_rotation(zen, az, phigeo, bfieldangle).T
    origin = numpy.array((cx, cy, cz))

    def GetXYZ(pos):
        return numpy.dot(Rt, pos) + origin
    return GetXYZ



######## oldschool
#def GetUVW(pos, cx, cy, cz, zen, az, phigeo, bfieldangle):
   #relpos = pos-numpy.array([cx,cy,cz])
   #inc=bfieldangle#/180.*numpy.pi #magnetic field direction on SKA site
   ##inc=-1.0554456843876574 #SKA
   ##B = numpy.array([0,numpy.cos(inc),-numpy.sin(inc)]) 
   #B = numpy.array([numpy.cos(phigeo)*numpy.sin(inc), numpy.sin(phigeo)*numpy.sin(inc),numpy.cos(inc)]) #from oliviers script including phigeo
   #v = numpy.array([numpy.cos(az)*numpy.sin(zen),numpy.sin(az)*numpy.sin(zen),numpy.cos(zen)])# *-1: account for angle conventions
   ##print v
   #vxB = numpy.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
   #vxB = vxB/numpy.linalg.norm(vxB)
   #vxvxB = numpy.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
   #return numpy.array([numpy.inner(v,relpos), numpy.inner(vxB,relpos),numpy.inner(vxvxB,relpos)]).T # vector dot