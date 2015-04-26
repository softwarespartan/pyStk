

import hfileLib;
import numpy as np;
import pyCoords;
import time;
import pyStk

def npvForStnList(hfile,stnList,type='glx'):
    
    N = 3*len(stnList);
    
    # mem alloc
    npv = np.arange(N)[:,np.newaxis]*np.nan;
    
    # extract the coords from hfile
    coords = hfile.coordsAsDict();
    
    # loop through the stnList and extract coords where possible
    for i in range(len(stnList)):
        
        # the i'th station name in the list
        stn = stnList[i].upper();
        
        # if the hfile contains this station then
        if coords.has_key(stn):
            
            # get coordinates from hfile in geocentric frame
            llr = np.array(coords[stn]['glx']);
            
            # convert from spherical to cartisian frame
            xyz = pyCoords.transform(llr,frm='spherical',to='cartisian')
            
            # that's it ... pack'em away
            npv[3*i:3*i+3] = xyz[:,np.newaxis];
                 
    return npv

def combine(npv1,npv2):
    
    tieIndx = np.nonzero(~np.isnan(npv1) & ~np.isnan(npv2))[0];
    
    npvSum = np.nansum(np.c_[npv1,npv2], 1);
    
    npvSum[tieIndx] = npvSum[tieIndx]/2.
    
    return npvSum[:,np.newaxis];

def stk():
    tstart  = time.time();
    
    hfiles = hfileLib.Hsrc(2010,103);
    hfiles.addProjectSrc('glbl','glx','../data/osuGLOBAL/gamit/orbit/solutions');
    hfiles.addProjectSrc('glbd','glx','../data/osuGLOBAL/gamit/densification/solutions');
    hfiles.addProjectSrc('capp','glx','../data/CAP/gamit/solutions');
    hfiles.addProjectSrc('anet','glx','../data/ANET/gamit/solutions');
    hfiles.addProjectSrc('gnet','glx','../data/GNET/gamit/solutions');
    
    # print the number of hfiles in the system
    print hfiles.size()
    
    # resolve any conflicting names
    hfiles.resolveNameConflicts();
    
    # get the stn list for npv syncronization
    stnNameList = list(set(hfiles.getStnNamesAsList()));
    
    npvF = npvForStnList(hfiles.hfileObjList[0],stnNameList);
    
    header    = 'file               status iter  RMS [mm]  wRMS [mm]  maxResid  %-out  N-out  N-ties  N-stns    Tx [mm]    Ty [mm]    Tz [mm]   Rx [nRad]   Ry[nRad]   Rz [nRad]'
    print header
    
    flist = list();
    
    for hfile in hfiles:
        hfileName = hfile.getUniqueHfileName();
    #    print 
    #    print
    #    print hfileName;
        
        # get the npv for the hfile
        npv = npvForStnList(hfile,stnNameList)
        
        # estimate transformation
        T,npvT,stats = pyStk.helmert(npv, npvF,percentile=30,reweight=5);
        
#        if np.any(np.abs(T['t'])>0.15) \
#            or np.any(np.abs(T['r'])>1e-7):
#            T,npvT,stats = pyStk.helmert(npv, npvGlobal);
            
        status = 'used'
        if stats['iter'] >= stats['maxIter']:
            status = 'failed'
            flist.append(hfileName+' failed to converge at '+str(stats['iter'])+' iterations')
        
        if stats['npts'] < 2:
            #print hfileName,'failed b/c only has ',stats['npts'],'common tie stations';
            #continue
            status = 'failed'
            flist.append(hfileName+' failed b/c only has '+str(stats['npts'])+' common tie stations')
            
        if np.any(np.abs(T['t'])>0.5) \
            or np.any(np.abs(T['r'])>1e-6):
            status = 'failed'
            flist.append(hfileName+' failed with crazy transform')
            
        
        # combine with new solution
        if status == 'used':
            npvF = combine(npvT.copy(),npvF.copy()).copy()
            
        if not np.any(np.abs(T['t'])>0.5) \
            and not np.any(np.abs(T['r'])>1e-6):
            T['t'] = T['t']/1e-3;
            T['r'] = T['r']/1e-9;
        else:
            T['t'] = np.ones_like(T['t'])*np.nan;
            T['r'] = np.ones_like(T['r'])*np.nan;
            
        if stats['RMS']/1e-3 > 100:
            stats['RMS'] = np.nan;
            
        if stats['wRMS']/1e-3 > 100:
            stats['wRMS'] = np.nan;
            
        id_status = '%s %6s %4d' % (hfileName,status,stats['iter']);
        rms       = '%8.2f %8.2f %11.2f' % (stats['RMS']/1e-3,stats['wRMS']/1e-3,stats['dvMax']/1e-3)
        pout      = '   %4.1f  %4d    %4d    %4d' % (stats['pout'],stats['nout'],stats['npts'],np.nonzero(~np.isnan(npv))[0].size/3);
        transform = '  %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f' % (T['t'][0],T['t'][1],T['t'][2],T['r'][0],T['r'][1],T['r'][2])
        print id_status,rms,pout,transform
        
    print
    for l in flist:
        print l
        
    print 'elapsed time: ', time.time() - tstart;
        
    return hfiles,stnNameList,npvF

stk()

