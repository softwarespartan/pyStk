

import numpy as np
import math

# global (?)
dr = math.pi/180.

class pyCoordsException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def __reshape(coords):
    
    #
    # input:  numpy array
    #
    # valid dimensions: 
    #     coord:  3 x 1
    #     npv:   3N x 1
    #     npvs    3 x N
    #
    # output:  numpy array of shape 3 x N
    #
    
    # first check that the number of elements is
    # divisible by 3.  if not, then not valid set
    if coords.size % 3 != 0:
        raise pyCoordsException('input coordinates must have 3N components');
    
    
    if coords.ndim == 1:
        coords = coords[:,np.newaxis]
        
    # figure out the actual dimensions of the coords
    n,m = coords.shape;
    
    # check for invalid 3N x M case
    if n != 3 and m !=1:
        raise pyCoordsException('input can not have shape 3N x M');
    
    # by default wed like to have 3 x N coord dimension
    # Only the npv form needs to be reshaped 
    if n != 3:
        coords = coords.reshape(3,n/3,order='F');
        
    return coords.copy();

def __restore(coords,n,m):
    
    # input 3 x N coords
    # output n x m coords
    
    # figure out the shape of the input
    dim_i,dim_j = coords.shape;
    
    # check for correct dimensions 
    if dim_i == n and dim_j == m:
        return coords;
    
    return coords.reshape(n,m,order='F')
    
def __x2e(x,y,z):
    
    # Converts global cartesian coordinates into  ellipsoidal (geodetic)
    # input:  x,y,z in meters
    # output: lat, lon, height in decimal degrees
    
    # 2x machine presision
    eps = 2.*10**-16;
    
    # convergence parameter
    tol=2.*eps;
    
    # max number of iterations
    itmax=10; 
    
    # Flattening in WGS-84
    f=1/298.257223563; 
    
    # WGS-84 equatorial radius in meters
    a=6378137.0;  
    
    esq=2.*f-f*f;
       
    p=math.sqrt(x**2+y**2); 
    
    # compute first approximation for lat
    lat=math.atan2(z/(1.-esq),p);
      
    #old lat   
    olat=lat;
    
    #iteration number
    it=0; 
    
    #inital difference
    diff=2.*tol;

    #start iterative process
    while diff > tol and it < itmax:
               
        N=a/math.sqrt(1. - esq*math.sin(lat)**2.);
        h=p/math.cos(lat) - N;                      
        lat=math.atan2(z,p*(1. - esq*N/(N+h)) );  
        
        diff=math.fabs(lat-olat); 
         
        olat=lat;
        it=it+1;
        
    #assign longitude    
    lon=math.atan2(y,x);    

    #check convergence
    if it >= itmax:
        raise pyCoordsException('transform from xyz to lat lon ht did not converge');

    # assign return list
    return lat*(1./dr),lon*(1./dr),h

def x2e(coords):
    
    # iterate over coords converting
    # from xyz to lat lon height
    
    # loop over the coords and transform
    for j in range(0,coords.shape[1]):
        coords[0,j],coords[1,j],coords[2,j] \
            = __x2e(coords[0,j],coords[1,j],coords[2,j])
            
    return coords;

def __e2x(lat,lon,ht):
    
    # convert ellipsoidal (geodetic) coordinates to global cartesian coordinates
    # input: lat,lon,height in decimal degrees
    # output: x,y,z in meters
    
    # convert to radians
    lat = lat*dr;
    lon = lon*dr;
    
    # Flattening according to WGS-84 
    f=1/298.257223563;  
    
    # Equatorial radius in meters. WGS-84 value (same as GRS 1980)   
    a=6378137.0;
    
    esq = 2.*f-f**2; 
     
    slat = math.sin(lat); clat = math.cos(lat); 
    slon = math.sin(lon); clon = math.cos(lon);
    
    # radius of curvature in prime vertical
    N = a/math.sqrt(1. - esq * slat**2);  
    
    # do it
    x = (N + ht) * clat * clon;
    y = (N + ht) * clat * slon;
    z = ( (1.- esq) * N + ht) * slat;
    
    # that's a [w]rap
    return x,y,z

def e2x(coords):
    
    # iterate over coords converting
    # from lat lon height to xyz 
    
    # loop over the coords and transform
    for j in range(0,coords.shape[1]):
        coords[0,j],coords[1,j],coords[2,j] \
            = __e2x(coords[0,j],coords[1,j],coords[2,j])

    return coords;

def __x2g(x,y,z):

    # transform xyz from cartesian coordinates to spherical coordinates lat,lon,r
    # 
    # input: x,y,z in meters 3 x M
    # output: lat,lon,ht in rad,rad,meters
    #
    
    lat = math.atan2(z,np.sqrt(x**2 + y**2));
    lon = math.atan2(y,x);
    ht  = math.sqrt(x**2 + y**2 + z**2);
    
    return lat*(1./dr),lon*(1./dr),ht;

def x2g(coords):
    
    # iterate over coords converting
    # from cartisian xyz to spherical lat,lon,r
    
    # loop over the coords and transform
    for j in range(0,coords.shape[1]):
        coords[0,j],coords[1,j],coords[2,j] \
            = __x2g(coords[0,j],coords[1,j],coords[2,j])
            
    return coords
            
def __g2x(lat,lon,r):

    # transform spherical lat,lon,r geographical coordinates 
    # to global cartesian xyz coordinates
    # 
    # input:  lat,lon,r in deg,deg,meters
    # output: x,y,z in meters 3 x M
    
    # convert to radians
    lat = lat * dr; lon = lon * dr;
    
    sla = math.sin(lat);  cla = math.cos(lat); 
    slo = math.sin(lon);  clo = math.cos(lon);

    x = r * cla * clo; 
    y = r * cla * slo; 
    z = r * sla;

    return x,y,z;

def g2x(coords):
    
    # iterate over coords converting
    # from spherical lat,lon,r to cartisian xyz
    
    # loop over the coords and transform
    for j in range(0,coords.shape[1]):
        coords[0,j],coords[1,j],coords[2,j] \
            = __g2x(coords[0,j],coords[1,j],coords[2,j])
            
    return coords
    
def __x2n(coords,ref):

    if ref.size != 3:
        raise pyCoordsException('reference position must be size 3 for __x2n');

    # convert the reference position to spherical
    ref = __x2g(ref[0],ref[1],ref[2]);
    
    sla=math.sin(ref[0]*dr); slo=math.sin(ref[1]*dr); 
    cla=math.cos(ref[0]*dr); clo=math.cos(ref[1]*dr);
    
    T = np.array([                             \
                    [-sla*clo, -sla*slo,  cla],\
                    [-slo,        clo,    0.0],\
                    [-cla*clo, -cla*slo, -sla] \
                 ]);
                  
    return T.dot(coords)
    
def x2n(coords,ref):

    # transform the coordinates of vectors V 
    # at reference points on the sphere
    # 
    # input:  coords 3 x M x,y,z in meters
    #            ref 3 x M or 3 x 1 x,y,z in meters
    #
    # output: coords 3 x M in local N, E, U coordinate system
    
    # control
    isSinglePoint = False;
    
    # make reference points into 3 x M;
    ref = __reshape(ref);

    if ref.size == 3:
        isSinglePoint = True;
    
    # at this point reference positions should match size
    if not isSinglePoint and ref.size != coords.size:
        raise pyCoordsException('reference position must be 3x1 or match size of input');

    if isSinglePoint:
        coords = __x2n(coords,ref);
    else:
        # loop over the coords and transform
        for j in range(0,coords.shape[1]):
            coords[0,j],coords[1,j],coords[2,j] \
                = __x2n(np.array([ coords[0,j], coords[1,j], coords[2,j] ]),\
                        np.array([    ref[0,j],    ref[1,j],    ref[2,j] ]));       
                        
    return coords;     
            
def __n2x(coords,ref):

    if ref.size != 3:
        raise pyCoordsException('reference position must be size 3 for __n2x');

    # convert the reference position to spherical
    ref = __x2g(ref[0],ref[1],ref[2]);
    
    sla=math.sin(ref[0]*dr); slo=math.sin(ref[1]*dr); 
    cla=math.cos(ref[0]*dr); clo=math.cos(ref[1]*dr);
    
    T = np.array([                             \
                    [-sla*clo, -sla*slo,  cla],\
                    [-slo,        clo,    0.0],\
                    [-cla*clo, -cla*slo, -sla] \
                 ]);
                  
    return T.transpose().dot(coords)
    
def n2x(coords,ref):

    # transform the coordinates of vectors V 
    # at reference points on the sphere
    # 
    # input:  coords 3 x M x,y,z in meters
    #            ref 3 x M or 3 x 1 x,y,z in meters
    #
    # output: coords 3 x M in local N, E, U coordinate system
    
    # control
    isSinglePoint = False;
    
    # make reference points into 3 x M;
    ref = __reshape(ref);

    if ref.size == 3:
        isSinglePoint = True;
    
    # at this point reference positions should match size
    if not isSinglePoint and ref.size != coords.size:
        raise pyCoordsException('reference position must be 3x1 or match size of input');

    if isSinglePoint:
        coords = __n2x(coords,ref);
    else:
        # loop over the coords and transform
        for j in range(0,coords.shape[1]):
            coords[0,j],coords[1,j],coords[2,j] \
                = __n2x(np.array([ coords[0,j], coords[1,j], coords[2,j] ]),\
                        np.array([    ref[0,j],    ref[1,j],    ref[2,j] ]));       
                        
    return coords;

def transform(coords,**kwargs):
    
    cartisian   = ('x','xyz','cart','cartisian');
    ellipsoidal = ('e','llh','ell','ellipsoidal');
    spherical   = ('g','llr','sph','spherical');
    local       = ('n','neu','local');
    
    validSrc = ('x','xyz','e','llh','n','neu','cart','cartisian','ell','ellipsoidal','local','g','llr','sph','spherical');
    validDst = ('x','xyz','e','llh','n','neu','cart','cartisian','ell','ellipsoidal','local','g','llr','sph','spherical');
    
    src = None; dst = None; ref = None;
    
    for key in kwargs:
        
        arg = kwargs[key];
        key = key.lower();
        
        if key == 'frm':
            if arg.lower() in validSrc:
                src = arg;
                
        elif key == 'to':
            if arg.lower() in validDst:
                dst = arg;
        
        elif key == 'at':
            ref = arg;
            
    if src == None or dst == None:
        raise pyCoordsException('must specify from and to for transform');
    
    if (src in local or dst in local) and ref == None:
        raise pyCoordsException('must provide reference postions when using local system')
            
    # check for proper input type
    if not isinstance(coords,np.ndarray):
        raise pyCoordsException('input must be numpy array')
    
    # figure out the original size\
    wasRowVector = False;
    if coords.ndim == 1:
        wasRowVector = True;
        coords = coords[:,np.newaxis];
    
    # record the shape of the original coords
    n,m = coords.shape;
    
    # convert coords to 3 x N shape
    try:
        coords = __reshape(coords);
    except pyCoordsException as e:
        raise e;
            
    # do the coordinate transformation
    if src in cartisian and dst in ellipsoidal:
        coords = x2e(coords);
        
    elif src in cartisian and dst in local:
        coords = x2n(coords,ref);
        
    elif src in cartisian and dst in spherical:
        coords = x2g(coords);   
        
    elif src in ellipsoidal and dst in cartisian:
        coords = e2x(coords);
        
    elif src in ellipsoidal and dst in local:
        coords = x2n(e2x(coords),ref)
        
    elif src in ellipsoidal and dst in spherical:
        coords = x2g(e2x(coords));
        
    elif src in local and dst in cartisian:
        coords = n2x(coords,ref);
        
    elif src in local and dst in ellipsoidal:
        coords = x2e(n2x(coords,ref));
        
    elif src in local and dst in spherical:
        coords = x2g(n2x(coords,ref));
    
    elif src in spherical and dst in cartisian:
        coords = g2x(coords);
        
    elif src in spherical and dst in ellipsoidal:
        coords = x2e(g2x(coords));
        
    elif src in spherical and dst in local:
        coords = x2n(g2x(coords),ref);
    
    # restore the shape of the original input
    if wasRowVector:
        return coords.flatten();
    else:
        return __restore(coords,n,m);
    