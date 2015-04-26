
import pyStk
import pyDate

#print pyStk.R(5, 6, 7)
#print
#print pyStk.dRdr1();
#print
#print pyStk.dRdr2();
#print
#print pyStk.dRdr3();

import scipy.io as sio
import numpy as np;
import time;

# create a time series with data
#ts = pyStk.pyTS().initFromMatFile('../data/ts.mat')

#npv1 = ts.npvForIndx(2005); 
#npv2 = ts.npvForIndx(2005);

#npv2[584,0]=0;

#T,npvT,stats = pyStk.helmert(npv1, npv2, percentile=30,reweightFactor=1.75,limit=2.5);
#
#print "final helmert transform params:"
#for k in T.keys():
#    print k,T[k].T
#
#print
#print 'iter:',stats['iter'];
#print
#print 'pout:',stats['pout'];
#print 'nout:',stats['nout'];
#print 'npts:',stats['npts'];
#print
#print ' RMS:',stats['RMS']/1e-3, '[mm]'
#print 'wRMS:',stats['wRMS']/1e-3,'[mm]';
#print
#print 'max resid:',stats['dvMax']/1e-3,'[mm]'
#print 'max resid indx:',stats['dvMaxIndx'][0]
#print 'max resid stn:',ts.stn_list[stats['dvMaxIndx']/3]
#
#print
#print 'ts size:',ts.size();
#print 'ts number of stations:',ts.numStns()
#print 'ts number of epochs',ts.numEpochs()

#dv = npv1 - npv2;
#weights = pyWeights.weights().computeWeights(dv[~np.isnan(dv)])


## init the stats
#iter = 0;
#rms  = 0;
#pout = 0;
#
## start the timer
#tstart = time.time();
#
## stk epoch per epoch
#for i in range(0,ts.size()-1):
#    
#    # get the data
#    npv1 = ts.npvForIndx(i);
#    npv2 = ts.npvForIndx(i+1);
#    
#    # stack the data 
#    T,npvT,stats = pyStk.helmert(npv1, npv2, percentile=20);
#    
#    # extract some statistics about the stk
#    iter += stats['iter'];
#    rms  += stats['wRMS'];
#    pout += stats['pout'];
#    
## stop the timer
#tstop = time.time();
#
## blab about it
#print 'total stk time: %4.1f seconds' % (tstop-tstart)
#print 'average iter: ', np.round(iter/ts.size());
#print 'average wRMS: %5.2f [mm]' % (rms/ts.size()/1e-3);
#print 'average pout: %5.1f' % (pout/ts.size())

tsFile = '/media/fugu/processing/projects/osuGLOBAL/shared/data/o11_anet_gnet_ts_NEW_PRIORS.mat';
tsFile = '/media/fugu/processing/projects/osuGLOBAL/shared/data/tmp.mat';
tsFile = '/Users/abelbrown/data/o11_napeos_itrf08_041212.mat'
ts = pyStk.pyTS().initFromMatFile(tsFile)
print ts.stn_list[5]

stnIndx =  ts.stnIndx(['ans_bir0','igs_yell'])
print stnIndx

for indx in stnIndx:
    print ts.stn_list[indx]

npvIndx = ts.npvIndx(['ans_bir0','igs_yell'])
print npvIndx

print 'does epoch 5 exists?',ts.epochExists(ts.epochs[5])
print 'what is the index of epoch 5?',ts.epochIndx(ts.epochs[5])
print 'what is the value of epoch 5?',ts.epochs[5]

print 'get npv for epoch 5'
print ts.npvForEpoch(ts.epochs[5]).shape;

print 'get npv for index 5'
print ts.npvForIndx(5).shape

print 'get xyz for igs_yell and igs_zimm'
print ts.xyzForStn(['igs_yell','igs_zimm']).shape

#for stn in ts.stns():
#    for xyz,d in ts.spvs(stn):
#        print stn,xyz,d.year,d.doy
