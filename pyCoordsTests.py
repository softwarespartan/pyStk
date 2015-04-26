

import pyStk;
import pyCoords;
import time
import numpy as np;
import math;

#ts = pyStk.pyTS().initFromMatFile('../data/ts.mat');
#npvT = ts.npvForIndx(2000);

#print 'before:', npvT.shape, npvT[576,0],npvT[577,0],npvT[578,0]
#print pyCoords.__x2e(npvT[576,0],npvT[577,0],npvT[578,0])
#
#npv = pyCoords.__reshape(npvT);
#print 'after: ',npv.shape, npv[0,192], npv[1,192], npv[2,192]
#print pyCoords.__x2e(npv[0,192], npv[1,192], npv[2,192])

#for npv,e in ts:
#    pyCoords.x2e(npv);
#    pyCoords.e2x(e);
    
#npv_x2g = pyCoords.x2g(npvT);
#print npv_x2g.shape
#print npvT[576,0],   npvT[577,0],   npvT[578,0]
#print npv_x2g[576,0],npv_x2g[577,0],npv_x2g[578,0]


#tstart = time.time();
#xyz = ts.xyzForStn(['igs_yell']);
#mean_xyz = np.mean(xyz,1);
#
#print xyz.shape;
#print mean_xyz.shape;
#npv_x2n = pyCoords.x2n(xyz, mean_xyz);
#tstop = time.time();
#print tstop - tstart


#print 'xyz: ', npvT[576,0],npvT[577,0],npvT[578,0];
#
#npv_e = pyCoords.transform(npvT, frm='cartisian', to='ellipsoidal')
#print 'x2e: ',npv_e[576,0],npv_e[577,0],npv_e[578,0]
#
#npv_x = pyCoords.transform(npv_e, frm='e', to='x')
#print 'e2x: ',npv_x[576,0],npv_x[577,0],npv_x[578,0]
#
#npv_n = pyCoords.transform(npv_x,frm='x',to='n',at=npvT);
#print 'x2n: ',npv_n[576,0],npv_n[577,0],npv_n[578,0]
#
#npv_n = pyCoords.transform(npv_e,frm='e',to='n',at=npvT);
#print 'e2n: ',npv_n[576,0],npv_n[577,0],npv_n[578,0]
#
#npv_nx = pyCoords.transform(npv_n,frm='n',to='x',at=npvT);
#print 'n2x: ',npv_nx[576,0],npv_nx[577,0],npv_nx[578,0]
#
#npv_ne = pyCoords.transform(npv_n,frm='n',to='e',at=npvT);
#print 'n2e: ',npv_ne[576,0],npv_ne[577,0],npv_ne[578,0]
#
#npv_ng = pyCoords.transform(npv_n,frm='n',to='g',at=npvT);
#print 'n2g: ',npv_ng[576,0],npv_ng[577,0],npv_ng[578,0]
#
#npv_gx = pyCoords.transform(npv_ng,frm='g',to='x',at=npvT);
#print 'g2x: ',npv_gx[576,0],npv_gx[577,0],npv_gx[578,0]

dr = math.pi/180;
lat = 0.7987232071688739  - float(0.5551942423196762E-08);
lon = -1.3626024290040994 + float(0.4905733872522222E-09);
r = 6367.3337427230617322 - float(0.1881525870011666E-04);
r = r/1e-3;
sph = np.array([lat*(1/dr),lon*(1/dr),r]);

print 'spherical coordinates: ',lat,lon,r
print pyCoords.transform(sph, frm='spherical', to='cartisian')


dr = math.pi/180;
lat = 0.8273480455582805  + float(0.2036669455549763E-07);
lon = -0.9194001780789284 - float(0.9838365675434765E-08);
r = 6366.6747206083946367 + float(0.3728246726314759E-04);
r = r/1e-3;
sph = np.array([lat*(1/dr),lon*(1/dr),r]);

print 'spherical coordinates: ',lat,lon,r
print pyCoords.transform(sph, frm='spherical', to='cartisian')


dr = math.pi/180;
lat = 0.8273480455582805     +float(0.1902355363654113E-07);
lon = -0.9194001780789284    -float(0.7743182595302539E-08);
r = 6366.6747206083946367     +float(0.2160005449332090E-04);
r = r/1e-3;
sph = np.array([lat*(1/dr),lon*(1/dr),r]);

print 'spherical coordinates: ',lat,lon,r
print pyCoords.transform(sph, frm='spherical', to='cartisian')

dr = math.pi/180;
lat = 0.8273480455582805     +float(0.1902355363654113E-07);
lon = -0.9194001780789284    -float(0.7743182595302539E-08);
r = 6366.6747206083946367     +float(0.2160005449332090E-04);
r = r/1e-3;
sph = np.array([lat*(1/dr),lon*(1/dr),r]);

print 'spherical coordinates: ',lat,lon,r
print pyCoords.transform(sph, frm='spherical', to='cartisian')