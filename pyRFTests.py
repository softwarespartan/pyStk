

import pyRF;
import pyStk;
import numpy as np;
import math

from scipy import stats;

#print rf.refData.shape
#print len(rf.refStnList)
#
#print rf.npv.shape;
#print rf.nvv.shape;
#print rf.refEpoch.shape;
#
#npv = rf.npvForEpoch(2003.50414524)

#print npv.shape

ts = pyStk.pyTS().initFromMatFile('../data/ts.mat');
rf = pyRF.pyRF('itrf08').initForStnList(map(str.lower,ts.stn_list));

fyear = ts.epochs[4000];

npv   = ts.npvForEpoch(fyear);
npvRF = rf.npvForEpoch(fyear);

print "mean: ", stats.nanmean(npv-npvRF)
print "median: ", stats.nanmedian(npv-npvRF)

print "Aligning epoch ",fyear
T,npvT,stats = pyStk.helmert(npv, npvRF);
print
print 'iter:',stats['iter'];
print
print 'pout:',stats['pout'];
print 'nout:',stats['nout'];
print 'npts:',stats['npts'];
print
print ' RMS:',stats['RMS']/1e-3, '[mm]'
print 'wRMS:',stats['wRMS']/1e-3,'[mm]';
print
print 'max resid:',stats['dvMax']/1e-3,'[mm]'
print 'max resid indx:',stats['dvMaxIndx'][0]
print 'max resid stn:',ts.stn_list[stats['dvMaxIndx']/3]


