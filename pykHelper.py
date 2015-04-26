

import hfileLib;
import numpy as np;
import pyCoords;
import time;
import pyStk
import sys,traceback;
import pyRF
import pyDate;
import stkData;
import pyWeights;

#np.seterr(all='raise');

def nanmean(npvs,weights=None,ref=None):
    
    if weights != None and (npvs.shape != weights.shape):
        print "npvs.shape =",npvs.shape, 'weights.shape = ',weights.shape
        raise pyStk.pyStkException('npvs and weights must be same shape for weighted mean');
    
    
    ###  BUG!!!  weights are init to one every time !!! #######
    
    # initialize the weights to same shape as npvs
    if weights == None:
        weights = np.ones(npvs.shape, dtype='float');
        
        # make sure that the nans in the weights match nans in npvs
        weights[np.isnan(npvs)]=np.nan;
    
    # construct the weighted mean
    npv = np.nansum(npvs*weights,1)/np.nansum(weights,1);
    
    if ref !=None:
        npv[np.nonzero(~np.isnan(ref))[0]] = ref[np.nonzero(~np.isnan(ref))[0]];
    
    # return as a column vector
    return npv[:,np.newaxis]

def findMergeMax(npvs):
    
    # number of possible columns
    N = npvs.shape[1];  M = N;
    
    # init the count/dist matrix to zero
    ncommon = np.zeros((N,M),dtype='float');
        
    # for each column find the number
    # of common stations with every other
    # colnmn including self.
    for i in range(0,N):
        npv_i = npvs[:,i];
        
        for j in range(0,M):
            
            count = 0;
            
            npv_j = npvs[:,j];
            
            if i >= j:
                ncommon[i,j]= np.nan;
                continue;
            
            # find non nan entries for both npvi npvsj
            # notice here we exploit npv construct
            indx = np.nonzero(~np.isnan(npv_i) & ~np.isnan(npv_j));
            
            if len(indx) != 0:
                count = indx[0].size/3.;
            else:
                count = 0;
                
            ncommon[i,j] = count;
            
    i,j = np.where(ncommon==np.nanmax(ncommon));
    
    #print ncommon
    
    return i[0],j[0],np.nanmax(ncommon);

def findMergeMaxForNpv(npvs,npv):
    
    N = npvs.shape[1];
    ncommon = np.zeros(N,dtype='float');
    npv = npv.flatten();
    
    for i in range(0,N):
        
        indx = np.nonzero(~np.isnan(npvs[:,i]) & ~np.isnan(npv))
            
        if len(indx) != 0:
                count = indx[0].size/3.;
        else:
            count = 0;
                
        ncommon[i] = count;
        
    j = np.where(ncommon==np.nanmax(ncommon));
    return np.array([j[0][0]]),np.nanmax(ncommon);

def nanMeanForRows(npvs,weights=None):
    
    if npvs == None:
        return None;
    
    if weights == None:
        weights = np.ones_like(npvs);
    
    # mem alloc
    rowMean = np.ones_like(npvs[:,0]);
    
    # calculate nanMedian for each row
    for i in range(0,np.shape(npvs)[0]):
        
        # get data for the i'th row
        rowData = npvs[i,:];
        
        # get weights for the i'th row
        rowWeights = weights[i,:];
        
        # figure out where the nan's are
        nanIndx = np.isnan(rowData);
        
        # compute nan median for this row data
        if np.all(nanIndx):
            rowMean[i] = np.nan;
        elif np.all(rowWeights[~nanIndx] == 0):
            rowMean[i] = np.nan;
        else:
            rowMean[i] = np.average(rowData[~nanIndx],weights=rowWeights[~nanIndx]);
            
    # make into column vector
    rowMean = rowMean[:,np.newaxis];         
    
    for j in range(0,np.shape(rowMean)[0]/3):
        
        if np.any(np.isnan(rowMean[3*j:3*j+3])):
            #print 'removing station at index =',j,rowMean[3*j:3*j+3].flatten() ;
            rowMean[3*j:3*j+3] = np.nan;
    
    return rowMean[:,np.newaxis];

def nanMedianForRows(npvs,incMin=np.Inf):
    
    if npvs == None:
        return None;
    
    # mem alloc
    rowMedian = np.ones_like(npvs[:,0]);
    
    # calculate nanMedian for each row
    for i in range(0,np.shape(npvs)[0]):
        
        # get data for the i'th row
        rowData = npvs[i,:];
        
        # figure out where the nan's are
        nanIndx = np.isnan(rowData);
    
        # compute the number of constituence
        inc = np.sum(~nanIndx);
        
        # compute nan median for this row data
        if np.all(nanIndx):
            rowMedian[i] = np.nan;
        elif inc < incMin:
            rowMedian[i] = np.nan;
        else:
            rowMedian[i] = np.median(rowData[~nanIndx]);
    
    return rowMedian;

class helper():
    
    def __init__(self,dataMgr):
        
        # translation and rotation tolerance
        self.ttol = 0.5
        self.rtol = 1e-6;
        
        # sigma tolerance 
        self.sigtol  = 1e-4;
        
        # maximum number of iterations
        self.maxIter = 20;
        
        # the object with initialized data sources 
        self.dataMgr       = dataMgr;
        
        # get the unique stn list for npv synchronization
        self.stnList = list(set(self.dataMgr.getStnNamesAsList()));
        
        # initialize each data src as npv for each column
        self.npvs          = self.mknpvs();
        
        # the initial target should be based on median (robust)
        # stations with single solution will not be in target
        self.npvTarget     = nanMedianForRows(self.npvs,2);
        
        # the init aligned npv's are just npvs
        self.npvsStacked   = self.npvs.copy();
        
        # everyone starts with weight = 1
        self.npvsWeights   = np.ones_like(self.npvs);
        
        # matrix of residuals for each npv w.r.t npvTarget
        self.npvsRedisuals = np.zeros_like(self.npvs);

        self.stats      = list();
        self.status     = list();
        self.transforms = list();
        self.summary    = list();                
        
        self.header    = 'Network                    status iter  maxResid  RMS [mm]  wRMS [mm]  IQR      %-out  N-out  N-ties  N-stns     Tx [mm]    Ty [mm]    Tz [mm]   Rx [nRad]   Ry[nRad]   Rz [nRad]  ZTD Bias[mm]'

        #self.summary.append(self.header);

    def __npv(self,dataObj):
        
        N = 3*len(self.stnList);
        
        # mem alloc
        npv = np.arange(N)[:,np.newaxis]*np.nan;
        
        # extract the coords from dataObj
        coords = dataObj.coordsAsDict();
        
        # loop through the stnList and extract coords where possible
        for i in range(len(self.stnList)):
            
            # the i'th station name in the list
            stn = self.stnList[i].upper();
            
            # if the dataObj contains this station then
            if coords.has_key(stn):
                
                # get coordinates from dataObj 
                xyz = np.array(coords[stn]);
                
                # that's it ... pack'em away
                npv[3*i:3*i+3] = xyz[:,np.newaxis];
                                     
        return npv;
    
    def __nav(self,dataObj):
        
        N = len(self.stnList);
        
        # mem alloc
        nav = np.arange(N)[:,np.newaxis]*np.nan;
        
        # extract the coords from dataObj
        atm = dataObj.atmAsDict();
        
        # loop through the stnList and extract coords where possible
        for i in range(len(self.stnList)):
            
            # the i'th station name in the list
            stn = self.stnList[i].upper();
            
            # if the dataObj contains this station then
            if atm.has_key(stn):
                
                # that's it ... pack'em away
                nav[i] = np.array(atm[stn]);
                     
        return nav
    
    def mknpvs(self):
        npvs = None;
        for dataObj in self.dataMgr:            
            npv = self.__npv(dataObj)
            if npvs == None:
                npvs = npv;
            else:
                npvs = np.c_[npvs,npv];
        return npvs;
    
    def __combine(self,npv1,npv2):
        
        tieIndx = np.nonzero(~np.isnan(npv1) & ~np.isnan(npv2))[0];
        
        npvSum = np.nansum(np.c_[npv1,npv2], 1);
        
        npvSum[tieIndx] = npvSum[tieIndx]/2.
        
        return npvSum[:,np.newaxis];
    
    def getStatus(self,T,stats):
        
        status = 'used'
        
        if stats['iter'] >= stats['maxIter']:
            status = 'failed'
        
        if stats['npts'] < 2:
            status = 'failed'
            
        if np.any(np.abs(T['t'])>self.ttol) \
            or np.any(np.abs(T['r'])>self.rtol):
            status = 'failed'
            
        if stats['pout'] > 60:
            status = 'failed';
            
        return status;
    
 
    def stk(self,**kwargs):
                        
        # get the data manager
        dataMgr     = self.dataMgr;
        
        # check if there is more than 1 network
        if len(dataMgr) == 1:
            return;
        
        iterCount   = 0;
        residualIQR = 0.025;
        
        # something very large
        previousResidualIQR = 1e10;
        
        # init with original
        previousNpvsStacked = self.npvsStacked.copy();
        
        # flag final iteration
        isFinalIter = False;
        
        # number of coordinates = numStns *3
        numCoords = np.sum(~np.isnan(self.npvs));
                
        while iterCount < self.maxIter:
            
            self.summary.append('iter = '+str(iterCount))
            self.summary.append(self.header);
            
            # align each column (npv) onto the npvTarget
            for j in range(0,self.npvs.shape[1]):
                                
                # align the j'th column onto the target
                try:
                    T,npvT,stats                            \
                        = pyStk.helmert(self.npvs[:,j],     \
                                        self.npvTarget,     \
                                        sigma=residualIQR   \
                          )
                      
                    # save the aligned npv    
                    self.npvsStacked[:,j] = npvT.flatten().copy();
                    
                    # save the weights associated with each network alignment
                    self.npvsWeights[:,j] = stats['weights'].flatten().copy();
                    
                except Exception, e:
                    #raise;
                    print dataMgr[j].name(),e;
                    self.npvsStacked[:,j] = np.ones_like(self.npvsStacked[:,j])*np.nan;
                    self.npvsWeights[:,j] = np.zeros_like(self.npvsStacked[:,j]);
                    
                # compute the residuals 
                self.npvsRedisuals[:,j] = self.npvTarget - self.npvsStacked[:,j];
                
                # figure out of alignment was successful 
                status = self.getStatus(T, stats);
                
                # if the network failed then zero it out via weights
                if status == 'failed':
                    self.npvsWeights[:,j] = np.zeros_like(self.npvsStacked[:,j]);
                    self.npvsWeights[np.isnan(self.npvs)] = np.nan;
                
                # record the stats and transformation parameters
                self.stats.append(stats);
                self.transforms.append(T);
            
                # record the summary for this network iteration
                self.summary.append(self.__mkSummary(dataMgr[j].name(), npvT, status, stats, T, np.NAN));
            
            
            # update the iteration count
            iterCount = iterCount + 1;
            
            # compute the residualIQR of all the residuals
            residualIQR = pyWeights.ipr(self.npvsRedisuals.flatten(), 50);
                        
            # compute difference between current and previous iteration IQRs
            diffIQR = previousResidualIQR - residualIQR;
            
            # compute the percentage of down weighted data
            pout = (np.sum(self.npvsWeights < 1)/(numCoords*1.0))*100.0;
            p1   = np.sum(self.npvsWeights[self.npvsWeights < 1 ] >= 0.5);
            p1   = (p1/(numCoords*1.0))*100;
            p2   = (np.sum(self.npvsWeights < 0.5)/(numCoords*1.0))*100;
            
            # generate some summary stats for this iteration
            self.summary.append( "pySTK: residualIQR = %6.2f [mm]" % (residualIQR/1e-3));
            self.summary.append( "pySTK: decrease in sigma = %6.2f [mm]" % (diffIQR/1e-3));
            self.summary.append( "pySTK: percentage of downweighted data = %4.1f, %4.1f, %4.1f" % (pout,p1,p2));
           
            # if this is the final iteration then we're done
            if isFinalIter:
                break;
            
            # check for convergence
            if  np.any((diffIQR) <self.sigtol):
                isFinalIter = True;  
                
            # special case that current iteration actually made IQR worse
            # if so roll back to previous iteration and boost the IQR sigma
            if diffIQR < 0 or pout > 35:
                
                # figure out how much to increase sigma
                increment = max(0.002, abs(diffIQR), (residualIQR+abs(diffIQR))/1.5);
                
                # restore previous iteration sigma plus fudge factor
                residualIQR      = previousResidualIQR + increment;
                
                # restore previous iteration aligned networks
                self.npvsStacked = previousNpvsStacked.copy();
            
                # blab about it
                self.summary.append( "pySTK: rolling back residualIQR = %6.2f [mm]" % (residualIQR/1e-3));
                
            # recompute the target from stacked npv
            self.npvTarget = nanMedianForRows(self.npvsStacked,2); 
            
            # save the stations for next round
            previousResidualIQR = residualIQR;
            previousNpvsStacked = self.npvsStacked.copy();
            
            # for readability
            self.summary.append("");

                
    def imposeRF(self,refFrame,year,doy,fid=None):
        
        # convert year doy to fractional year
        fyear = pyDate.yeardoy2fyear(year, doy);
        
        # need to translate the 4 char names of 
        #current npvs to real stnIDs 8 char
        stnIdList = list();
        translate = self.dataMgr.getReverseAliasMap();
        for i in range(len(self.stnList)):
            stnId = list(translate[self.stnList[i]])[0];
            stnId = stnId.lower().replace('::','_');
            stnIdList.append(stnId);
        
        # initialize new reference frame
        rf = pyRF.pyRF(refFrame).initForStnList(stnIdList);
        
        npvTarget = rf.npvForEpoch(fyear);
        
        T,npvT,stats = pyStk.helmert(self.npvF, npvTarget,reweightFactor=5,limit=1.75,maxIter=100);
        
        status = self.getStatus(T, stats);
        
        summary = self.__mkSummary(refFrame, npvT, status, stats, T, np.NAN);
        
        if not fid == None:
            fid.write('\n');
            fid.write('Reference Frame:\n');
            fid.write(summary+'\n')
        else:
            print
            print 'Reference Frame:'
            print summary;
            
        if status == 'used':
            self.npvF = npvT;
                        
    def __mkSummary(self,hfileName,npv,status,stats,T,atmBias):
                        
        if not np.any(np.abs(T['t'])>self.ttol*10) \
            and not np.any(np.abs(T['r'])>self.rtol*10):
            T['t'] = T['t']/1e-3;
            T['r'] = T['r']/1e-9;
        else:
            T['t'] = np.ones_like(T['t'])*np.nan;
            T['r'] = np.ones_like(T['r'])*np.nan;
            
        if stats['RMS']/1e-3 > 1000:
            stats['RMS'] = np.nan;
            
        if stats['wRMS']/1e-3 > 1000:
            stats['wRMS'] = np.nan;
            
        if stats['dvMax']/1e-3 > 1000:
            stats['dvMax'] = np.nan;
            
        iqr = stats['iqr']/1e-3;
            
        id_status = '%-25s %6s %4d' % (hfileName,status,stats['iter']);
        rms       = ' %8.2f %8.2f %8.2f %9.2f' % (stats['dvMax']/1e-3,stats['RMS']/1e-3,stats['wRMS']/1e-3,iqr)
        pout      = '    %5.1f  %4d    %4d    %4d' % (stats['pout'],stats['nout'],stats['npts'],np.nonzero(~np.isnan(npv))[0].size/3);
        transform = '  %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f' % (T['t'][0],T['t'][1],T['t'][2],T['r'][0],T['r'][1],T['r'][2]);
        atm       = '    %7.3f' % (atmBias/1e-3);
        
        return id_status+" "+rms+" "+pout+" "+transform+" "+atm;           
            
    def printSummary(self,fid=None):
        
        if not fid == None:
            for s in self.summary:
                fid.write(s+'\n')
        else:
            for s in self.summary:
                print s;
            
    def isValid(self):
        
        if self.npvs != None:
            return True;
        else:
            return False;
            
    def printCoords(self,fid=None):
        
        header = '   stn     #ref         X                  Y                  Z           --  std before [mm]  --     --  std after [mm]  --    ZTD [m]     -- std [mm] --';
        translate = self.dataMgr.getReverseAliasMap();

        #npvsF = nanmean(self.npvsStacked,None);
        npvsF = nanMeanForRows(self.npvsStacked,self.npvsWeights);
        
        if fid != None:
            fid.write('\n');
            fid.write(header+'\n');
        else:
            print
            print header;

        for i in range(len(self.stnList)):
            
            stnId = list(translate[self.stnList[i]])[0];
            xyz   = npvsF[3*i:3*i+3,:];
            #atm   = self.navF[i];
            atm   = np.NAN;
            
            # compute std before alignment 
            xyzStd = np.zeros(3);
            refcnt = 0;
            for j in range(3):
                
                xi = self.npvs[3*i+j,:].flatten();
                indx = np.nonzero(~np.isnan(xi))[0];
                refcnt = indx.size;
                if indx.size <= 1:
                    continue
                xyzStd[j] = np.std(xi[indx])/1e-3;
                
            # compute std after alignment    
            xyzStdT = np.zeros(3);
            for j in range(3):
                xi = self.npvsStacked[3*i+j,:].flatten();
                indx = np.nonzero(~np.isnan(xi))[0];
                refcnt = indx.size;
                if indx.size <= 1:
                    continue
                xyzStdT[j] = np.std(xi[indx])/1e-3;
                
            ################# ATM #######################
            
#            # compute std ATM before and after adjustment 
#            atmStd  = 0.0;
#            atmStdT = 0.0;
#            atmRefCount = 0;
#            
#            # get the atm data before and after adjustment
#            ai = self.navs[i,:].flatten();
#            aiT = self.navsF[i,:].flatten();
#            
#            # figure out where the data is
#            indx = np.nonzero(~np.isnan(xi))[0];
#            
#            # compute the number of solutions that 
#            # reference this station
#            atmRefCount = indx.size;
#            
#            # make sure there exists a reference 
#            if atmRefCount > 1:
#                # finally comput the stand dev.
#                atmStd  = np.std(ai[indx])/1e-3;
#                atmStdT = np.std(aiT[indx])/1e-3;
            
            atm = np.NAN; atmStd = np.NAN; atmStdT = np.NAN;
            
            # make the string for printing to stdout or file
            coordStr = '%9s %3d %18.8f %18.8f %18.8f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %10.5f %8.2f %8.2f' % (stnId,refcnt,xyz[0,:],xyz[1,:],xyz[2,:],xyzStd[0],xyzStd[1],xyzStd[2],xyzStdT[0],xyzStdT[1],xyzStdT[2],atm,atmStd,atmStdT)

            # print the string;
            if fid != None:
                fid.write(coordStr+'\n');
            else:
                print coordStr
                
    def export(self, path='.'):
        
        # init
        ts = dict();
        
        # export with full station ids
        translate = self.dataMgr.getReverseAliasMap();
        
        # init the station list
        ts['stnm'] = list();
        
        for stnName in self.stnList:
            stnId = list(translate[stnName])[0];
            ts['stnm'].append(stnId);
            
           
        # compute the fractional year 
        ts['epochs']   = pyDate.Date(year=self.dataMgr.year,doy=self.dataMgr.doy).fyear;
        
        # the actual data
        ts['npvs'] = self.npvs;
        
        # the pyk daily merge
        #npvsF = nanMeanForRows(self.npvsStacked,self.npvsWeights);
        #ts['npv_pyk'] = npvsF;
        
        # generate a list of network names
        net_list= list();
        for obj in self.dataMgr:
            net_list.append(obj.name());
        ts['net_list'] = net_list;
        
        # import mat file capabilities
        import scipy.io;
        import os;
        
        # make sure the path given actually exists
        if not os.path.isdir(path):
            raise pyStk.pyStkException('path '+path+' does not exist\n');
        
        # generate the tsd file name
        fileName = 'tsd_'+str(self.dataMgr.year)+'_'+str(self.dataMgr.doy)+'.mat';
        
        # create the full export file path with file name
        filePath = os.path.join(path,fileName);
        
        # do it 
        scipy.io.savemat(filePath, ts, oned_as='column', do_compression=True)
            
                
if __name__ == '__main__':
    for year in range(2012, 2013):
        for doy in range(1,366):
        
            #year = 2010;
            #doy  = 100;
            
            try:
            
                dataMgr = stkData.DataMgr(year,doy);
                
                #snxData = stkData.SnxSrc(year,doy);
                #snxData.addProjectSrc('esa2','es2','snx','/media/fugu/processing/products/esa_repro2')
                #snxData.addProjectSrc('glbl','o20','snx','/media/fugu/processing/projects/osuGLOBAL/napeos/orbit/solutions');
                #snxData.addProjectSrc('glbd','o19','snx','/media/fugu/processing/projects/osuGLOBAL/napeos/densification/solutions');
                #snxData.addProjectSrc('glbd','o17','snx','/media/fugu/processing/projects/osuGLOBAL/napeos/densification/solutions');
                #snxData.addProjectSrc('glbd','osf','snx','/media/fugu/processing/projects/osuGLOBAL/napeos/orbit/solutions');
                #dataMgr.addData(snxData);
                
                #hfileData = stkData.HfileSrc(year,doy);
                #hfileData.addProjectSrc('glbl','gfx','glx','/media/fugu/processing/projects/osuGLOBAL/gamit/orbit/solutions/');
                #hfileData.addProjectSrc('glbl','os5','glx','/media/fugu/processing/projects/osuGLOBAL/gamit/orbit/solutions/');
                
                #hfileData.addProjectSrc('glbl','osf','glx','/media/fugu/processing/projects/osuGLOBAL/gamit/orbit/solutions/');
                #hfileData.addProjectSrc('glbd','osf','glx','/media/fugu/processing/projects/osuGLOBAL/gamit/densification/solutions/');
                #dataMgr.addData(hfileData);
                
                prtData = stkData.PrtSrc(year,doy);
                prtData.addProjectSrc('glbk', 'gk', 'prt', '/media/fugu/processing/projects/PRODUCT/globk/solutions/');
                dataMgr.addData(prtData);
                
                
                if len(dataMgr) == 0:
                    print year,doy, 'has no data'
                    continue;
                
                # resolve any conflicting names
                dataMgr.resolveNameConflicts();
                    
                hm = helper(dataMgr);
                #fh = hm.stk();
                #hm.printSummary()
                #hm.imposeRF('itrf08', year, doy)
                #hm.printCoords();
                
                #fid.close();
                print 'exporting ',year,doy;
                #hm.export('/Users/abelbrown/data/osf/gamit/orbit_densification');
                hm.export('/media/fugu/processing/projects/PRODUCT/globk/mat');
                
            except Exception, e:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                
                print "*** print_tb:"
                traceback.print_tb(exc_traceback, limit=10, file=sys.stdout)
                print
                print "*** print_exception:"
                traceback.print_exception(exc_type, exc_value, exc_traceback,limit=20, file=sys.stdout)
            
            
                
    
