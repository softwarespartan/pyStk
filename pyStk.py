'''
Created on May 30, 2011

@author: abel
'''
import numpy as np
import numpy.linalg;
import pyWeights;
import pyDate;

class pyTSException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class pyStkException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
                    
__dRdr1=np.array(                         \
                      [                   \
                        [  0,  0,  0 ],   \
                        [  0,  0,  1 ],   \
                        [  0, -1,  0 ]    \
                      ]                   \
                    );
                 
__dRdr2 = np.array(                      \
                      [                  \
                        [  0,  0, -1 ],  \
                        [  0,  0,  0 ],  \
                        [  1,  0,  0 ]   \
                      ]                  \
                    );
  
__dRdr3 = np.array(                      \
                      [                  \
                        [  0,  1,  0 ],  \
                        [ -1,  0,  0 ],  \
                        [  0,  0,  0 ]   \
                      ]                  \
                    );

def R(r):
    return np.array(                               \
                      [                            \
                        [    1,   r[2],  -r[1] ],  \
                        [ -r[2],    1,    r[0] ],  \
                        [  r[1], -r[0],     1  ]   \
                      ]                            \
                    );
                    
def applyRotation(npv,r):
    
    if r.size != 3:
        raise pyStkException('rotation input array must contain exactly 3 values');
    
    # change from 3N x 1 to 3 x N vector
    npv = npv.reshape(npv.size/3,3).transpose()
    
    # rotate the coordinates
    rnpv =  R(r).dot(npv);
    
    # convert back to column vector and return
    return rnpv.transpose().reshape(rnpv.size,1);

def applyTranslation(npv,t):
    
    if t.size != 3:
        raise pyStkException('translation input array must contain exactly 3 values');
    
    # make sure t  is a column vector
    if t.shape == (3,):
        t = t[:, np.newaxis];
        
    # tile the 3 x 1 vector to match size of npv
    t = np.tile(t,(npv.size/3,1));
    
    return npv - t;
    
def applyScale(npv,s):
    
    if s.size != 1:
        raise pyStkException('scale input must contain/be exactly 1 value');
    
    return npv * ( 1 + s );
    
def applyTransform(npv, **kwargs):
        
    for key in kwargs:
        
        if key.lower() == 'transform':
            
            params = kwargs[key];
            
            if isinstance(params,dict):                
                if params.has_key('r'):
                    npv = applyTransform(npv, rotation=params['r']);
                if params.has_key('t'):
                    npv = applyTransform(npv, translation=params['t']);
                if params.has_key('s'):
                    npv = applyTransform(npv, scale=params['s']);
                        
            elif isinstance(params,numpy.ndarray):
                if params.size >=3:
                    npv = applyTransform(npv,rotation=params[0:3])
                if params.size >= 6:
                    npv = applyTransform(npv,translation=params[3:6]);
                if params.size == 7:
                    npv = applyTransform(npv,scale=params[6]);
            else:
                raise pyStkException('transform type must be numpy.ndarray or dict')
        
        if key.lower() in ('t','trans','translation'):
            npv = applyTranslation(npv,kwargs[key]);
            
        elif key.lower() in ('r','rot','rotation'):
            npv = applyRotation(npv,kwargs[key]);
            
        elif key.lower() in ('s','scale'):
            npv = applyScale(npv,kwargs[key]);
            
    return npv

def jvi(v):
    
    # the simplified readable version of the 
    # jacobian for i'th vector v in the npv

    x = v[0];
    y = v[1];
    z = v[2];
    
    return np.array(                                 \
                    [                                \
                     [ 0., -z,  y, -1.,  0.,  0., x],\
                     [ z,  0., -x,  0., -1.,  0., y],\
                     [-y,  x,  0.,  0.,  0., -1., z] \
                    ]                                \
                   )

def jacobian(npv):

    # compute the number of coordinates
    N = npv.size;
    
    # reshape the npv 3N x 1 into matrix 3 x N
    npvMat = npv.reshape(N/3,3).transpose(); 
    
    # multiply wrt derivative of rotation matrix
    # and then turn into a column vectors
    c1 = __dRdr1.dot(npvMat).transpose().reshape(N,1);
    c2 = __dRdr2.dot(npvMat).transpose().reshape(N,1);
    c3 = __dRdr3.dot(npvMat).transpose().reshape(N,1);
    
    return np.hstack(                              \
                     (                             \
                      c1, c2, c3,                  \
                      np.tile(-np.eye(3),(N/3,1)), \
                      npv                          \
                     )                             \
                    );

def helmert(npv,npvTarget,**kwargs):
    
    # convert to column vectors
    npv       = npv[:,np.newaxis];
    npvTarget = npvTarget[:,np.newaxis];
    
    # first things first need to check shapes
    if npv.shape !=npvTarget.shape:
        raise pyStkException('npv and target must be same shape');
    
    # next init the jacobian and initial 
    # parameter guess for the transformation
    # note: here assume not estimating scale
    J      = jacobian(npv)[:,0:6];
    params = np.array([0., 0., 0., 0., 0., 0.])[:, np.newaxis];
    
    # make copy of npv for final npvT computation
    # so that the output npvT matches input npv with nans
    npvOri = np.copy(npv); npvTargetOri = np.copy(npvTarget);
    
    # init weights
    weights = np.ones_like(npv);
    
    #max number of iterations
    maxIter = 30;
    
    #control variables (assume meters for units)
    iter                  = 1;
    converged             = False;
    rotation_tolerance    = 5e-11;
    translation_tolerance = 1e-4;
    wRMS_tolerance        = 1e-3;
    shouldComputeWeights  = True;
    shouldComputeStats    = True;
    
    # control variables for residual reweighting
    percentile     = 50;
    limit          = 1.75;
    reweightFactor = 4;
    aprioriSigma   = None;

    # check for key word args to over ride default settings
    for key in kwargs:

        if key.lower() in ('withscale','usescale','estscale','estimatescale'):
            if str(kwargs[key]).lower() in ('yes','true','1'):
                J = jacobian(npv);
                params = np.array([0., 0., 0., 0., 0., 0., 0.])[:, np.newaxis];
        
        if key.lower() in ('w','weights'):
            weights = kwargs[key];
            shouldComputeWeights = False;
            if weights.shape != npv.shape:
                raise pyStkException('weights and npv must be same shape');
            
        if key.lower() in ('maxiter','maxiteration','itermax', 'iterationmax'):
            maxIter = kwargs[key];
            if maxIter < 1:
                raise pyStkException('number of iterations must be at least 1');
                
        if key.lower() in ('rot_tol','rotation_tolerance','rtol'):
            rotation_tolerance = kwargs[key];
            if rotation_tolerance < 0:
                raise pyStkException('rotation tolerance must be positive');
            
        if key.lower() in ('trans_tol','translation_tolerance','ttol'):
            translation_tolerance = kwargs[key];
            if translation_tolerance < 0:
                raise pyStkException('translation tolerance must be positive');
            
        if key.lower() == 'percentile':
            percentile = kwargs[key];
            if percentile <=0:
                raise pyStkException('weight percentile must be > 0');
        
        if key.lower() == 'limit':
            limit = kwargs[key];
            if limit <= 0:
                raise pyStkException('weight overage limit must be > 0');
            
        if key.lower() == 'reweightfactor':
            reweightFactor = kwargs[key];
            if reweightFactor <=0:
                raise pyStkException('reweightFactor must be > 0');
            
        if key.lower() in ('computestats','shouldcomputestats','stats'):
            if str(kwargs[key]).lower() in ('no','false','0'):
                shouldComputeStats = False;
                
        if key.lower() in ('aprioriSigma','apriorisigma','sigma'):
            aprioriSigma = kwargs[key];
            if aprioriSigma != None and aprioriSigma <= 0:
                raise pyStkException('a priori aprioriSigma must be greater than zero');
   
    # lastly, need to remove any nan values from the calculation
    indx      = np.nonzero( ~np.isnan(npv) & ~np.isnan(npvTarget) )[0];
    npv       =       npv[indx];
    npvTarget = npvTarget[indx];
    weights   =   weights[indx];
    J         =         J[indx,:];
    
    # make sure there are some points left to align with
    if npv.size/3 == 0:
        raise pyStkException('pystk: NPTS = 0, no common points to align with')

    # init transformed npv, residuals and weights
    npvT = npv;
    dv   = npvTarget - npvT;
    
    # initialize the weights generator
    weightsFactory            = pyWeights.weights();
    weightsFactory.percentile = percentile;
    weightsFactory.limit      = limit;
    weightsFactory.reweight   = reweightFactor;
    
    # set a priori sigma if defined
    if aprioriSigma != None:
        weightsFactory.aprioriSigma = aprioriSigma;
    
    # init weights if needed
    if shouldComputeWeights:
        weights = weightsFactory.computeWeights(dv);
        
    # iterate for convergence         
    while iter <= maxIter and not converged:
        
        # update iteration count
        iter = iter + 1;
        
        # finally(!) some inversion ... 
        try:  
            drts = numpy.linalg.lstsq(J*weights,dv*weights)[0];
        except:
            # pass it along, nothing we can do heres
            raise 
            
        # check if i'th set of rotations and translations 
        # estimates below tolerance i.e. check for convergence
        if iter > 1 and np.any( np.abs( drts[0:3]) < rotation_tolerance) \
                and np.any( np.abs( drts[3:6]) < translation_tolerance):
            converged = True;
            break;
        
        # assign update
        params = params + drts;
        
        # apply update
        npvT = applyTransform(npv,transform=params);
        
        # compute the difference between transform and target
        dv = npvTarget - npvT;
        
        # compute/recompute weights
        if shouldComputeWeights:
            weights = weightsFactory.computeWeights(dv);
          
        # recompute the weighted RMS for the residuals  
        wRMS = np.linalg.norm(dv*weights)/np.sqrt(dv.size);
        
#        # if wRMS is very small additional iterations could 
#        # cause strange floating point and overflow exceptions
#        if wRMS < wRMS_tolerance:
#            converged = True;
#            break;
        
    # reconstruct npvT in original size
    npvT = applyTransform(npvOri,transform=params)

    # pkg up statistics and config stuff
    stats = dict();
     
    if shouldComputeStats:
        
        stats['RMS']       = np.linalg.norm(dv)/np.sqrt(dv.size);
        stats['wRMS']      = np.linalg.norm(dv*weights)/np.sqrt(dv.size);
        
        stats['npts']      = npv.size/3;
        stats['nout']      = np.nonzero(weights<1)[0].size;
        
        stats['weights']   = weightsFactory.computeWeights(npvTargetOri-npvT);
        stats['pout']      = (np.nonzero(weights<1)[0].size/float(weights.size))*100;
        
        stats['iqr']       = pyWeights.ipr(dv, 50);
        stats['dvMax']     = np.nanmax(np.abs(dv));
        stats['dvMaxIndx'] = np.nonzero(np.abs(npvTargetOri-npvT)==stats['dvMax'])[0]
            
        stats['iter']                  = iter;
        stats['maxIter']               = maxIter;
        stats['converged']             = converged;
        stats['rotation_tolerance']    = rotation_tolerance;
        stats['translation_tolerance'] = translation_tolerance;
        
        stats['limit']          = limit;
        stats['percentile']     = percentile;
        stats['reweightFactor'] = reweightFactor;
    
    T      = dict();
    T['r'] = params[0:3];
    T['t'] = params[3:6];
    
    if params.size == 7:
        T['s'] = params[6];
    
    return T,npvT,stats

def estAtmBias(atm,atmTarget,**kwargs):
    
    # NOTE:  atm in meters
    
    # first things first need to check shapes
    if atm.shape !=atmTarget.shape:
        raise pyStkException('atm and target must be same shape');
    
    # next init the jacobian
    # notice that J*v = 1*atm = atm
    J      = atm;
    atmAdj = np.array([0.0]);
    
    # make copy of atm for final atmT computation
    # so that the output atmT matches input atm with nans
    atmOri = np.copy(atm); atmTargetOri = np.copy(atmTarget);
    
    # init weights
    weights = np.ones_like(atm);
    
    #max number of iterations
    maxIter = 20;
    
    #control variables (assume meters for units)
    iter = 1;
    converged = False;
    tolerance = 1e-5;
    shouldComputeWeights = True;
    shouldComputeStats   = True;
    
    # control variables for residual reweighting
    percentile     = 50;
    limit          = 1.65;
    reweightFactor = 1.5;

    # check for key word args to over ride default settings
    for key in kwargs:
        
        if key.lower() in ('w','weights'):
            weights = kwargs[key];
            shouldComputeWeights = False;
            if weights.shape != atm.shape:
                raise pyStkException('weights and atm must be same shape');
            
        if key.lower() in ('maxIter','maxIteration','iterMax', 'iterationMax'):
            maxIter = kwargs[key];
            if maxIter < 1:
                raise pyStkException('number of iterations must be at least 1');
                
        if key.lower() in ('tol','tolerance'):
            tolerance = kwargs[key];
            if tolerance < 0:
                raise pyStkException('tolerance must be positive');
            
        if key.lower() == 'percentile':
            percentile = kwargs[key];
            if percentile <=0:
                raise pyStkException('weight percentile must be > 0');
        
        if key.lower() == 'limit':
            limit = kwargs[key];
            if limit <= 0:
                raise pyStkException('weight overage limit must be > 0');
            
        if key.lower() == 'reweightfactor':
            reweightFactor = kwargs[key];
            if reweightFactor <=0:
                raise pyStkException('reweightFactor must be > 0');
            
        if key.lower() in ('computestats','shouldcomputestats','stats'):
            if str(kwargs[key]).lower() in ('no','false','0'):
                shouldComputeStats = False;
   
    # lastly, need to remove any nan values from the calculation
    indx      = np.nonzero( ~np.isnan(atm) & ~np.isnan(atmTarget) )[0];
    atm       = atm[indx];
    atmTarget = atmTarget[indx];
    weights   = weights[indx];
    J         = J[indx,:];
    
    if atm.size == 0:
        raise pyStkException('NPTS = 0, no common points to est bias')

    # init transformed atm
    atmT = atm;
    
    # initialize the weights generator
    weightsFactory            = pyWeights.weights();
    weightsFactory.percentile = percentile;
    weightsFactory.limit      = limit;
    weightsFactory.reweight   = reweightFactor;
        
    # iterate for convergence         
    while iter <= maxIter and not converged:
        
        # compute the difference between atm and target atm
        dv = atmTarget - atmT;
        
        # compute/recompute weights
        if shouldComputeWeights :
            weights = weightsFactory.computeWeights(dv)
        
        # finally(!) some inversion ... 
        try:  
            adj = numpy.linalg.lstsq(J*weights,dv*weights)[0];
        except:
            # pass it along, nothing we can do heres
            raise 
            
        # check if i'th set of rotations and translations 
        # estimates below tolerance i.e. check for convergence
        if iter > 1 and adj <= tolerance :
            converged = True;
            break;
        
        # assign update
        atmAdj = atmAdj + adj;
        
        # apply update
        atmT = atm + atmAdj;
        
        # update iteration count
        iter = iter + 1;
        
    # reconstruct atmT in original size
    atmT = atmOri + atmAdj;

    # pkg up statistics and config stuff
    stats = dict();
     
    if shouldComputeStats:
        
        stats['RMS']       = np.linalg.norm(dv)/np.sqrt(dv.size);
        stats['wRMS']      = np.linalg.norm(dv*weights)/np.sqrt(dv.size);
        
        stats['npts']      = atm.size;
        stats['nout']      = np.nonzero(weights<1)[0].size;
        
        stats['weights']   = pyWeights.weights().computeWeights(atmTargetOri-atmT);
        stats['pout']      = (np.nonzero(weights<1)[0].size/float(weights.size))*100;
        
        stats['dvMax']     = np.nanmax(np.abs(dv));
        stats['dvMaxIndx'] = np.nonzero(np.abs(atmTargetOri-atmT)==stats['dvMax'])[0]
            
        stats['iter']      = iter;
        stats['maxIter']   = maxIter;
        stats['converged'] = converged;
        stats['tolerance'] = tolerance;

        stats['limit']          = limit;
        stats['percentile']     = percentile;
        stats['reweightFactor'] = reweightFactor;
    
    T = dict(); T['t'] = atmAdj;
    
    return T,atmT,stats

class pyTS():
    
    def __init__(self):
    
        self.npvs     = np.array([]);
        self.epochs   = np.array([]);
        self.stn_list = list();
        
        self.iterIndx = None;

    def initFromMatFile(self,matFile):
        
        import os;
        import scipy.io;
        
        if not os.path.isfile(matFile):
            raise pyTSException('matFile '+matFile+' not found');
        
        try:
            # open as hd5 file mat file v7.3
            data = scipy.io.loadmat(matFile);
            
            print data.keys();
            
            # snag the data
            if 'stnList' in data.keys():
                self.stn_list = data['stnList'];
            else:
                self.stn_list = data['stn_list'];
                
            self.epochs   = data['epochs'];
            self.npvs     = data['npvs'];
            
            if 'nsvs' in data.keys():
                self.npvs_sigma = data['nsvs'];
            elif 'npvs_sigmas' in data.keys():
                self.npvs_sigma = data['npvs_sigmas'];
            elif 'npvs_sigma' in data.keys():
                self.npvs_sigma = data['npvs_sigma'];
            elif 'npvsSigma' in data.keys():
                self.npvs_sigma = data['npvsSigma'];
            elif 'npvsSigmas' in data.keys():
                self.npvs_sigma = data['npvsSigmas'];
            else:
                self.npvs_sigma = None;
            
        except Exception as e:
            os.sys.stderr.write(str(e)+'\n');
            raise pyTSException('error reading mat file as ts structure');
        
        # clean up the data a bit
        self.epochs.resize(self.epochs.size)
        self.stn_list = [str(s[0][0]) for s in self.stn_list];
        
        return self
    
    def size(self):
        return self.npvs.shape[1];
    
    def numStns(self):
        return len(self.stn_list);
    
    def numEpochs(self):
        return self.epochs.size;
        
    def epochExists(self,epoch):
           
        indx = self.epochIndx(epoch);
                    
        if indx.size != 0:
            return True;
        else:
            return False;
        
    def stnExists(self,stn):
        try:
            self.stn_list.index(stn.upper());
            return True;
        except ValueError:
            return False;    
        
    def stnIndx(self,stns):
        #
        # returns the stn_list index for a given station key(s) as numpy array
        #
        # NOTE: stations that do not exist are ignored
        #
        
        # making sure you assholes get the input right
        if not isinstance(stns,str) and not isinstance(stns,list):
            raise pyTSException('unrecognized data type as input');
        
        # turn into list so we can iterate over it
        if isinstance(stns,str):
            stns = [stns];
        
        # init index array (as type integer for indexing)
        indx = np.array([],dtype=int);
        
        # loop over stations
        for s in stns:
            s = s.upper();
            try:
                indx = np.r_[indx,self.stn_list.index(s)]
            except ValueError:
                indx = indx          

        return indx
    
    def npvIndx(self,stns):

        # making sure you muther fuckers get the input right
        if not isinstance(stns,str) and not isinstance(stns,list):
            raise pyTSException('unrecognized data type as input');
        
        # turn into list so we can iterate over it
        if isinstance(stns,str):
            stns = [stns];
        
        # compute the stn list indices for each station
        stnIndx = self.stnIndx(stns);
        
        # mem alloc
        indx = np.zeros(3*stnIndx.size,dtype=int);
        
        # loop over and set indicies explicitly
        for i in range(0,stnIndx.size):
            indx[3*i:3*(i+1)] = np.array([3*stnIndx[i],3*stnIndx[i]+1,3*stnIndx[i]+2]); 
        
        # simple as that
        return indx
    
    def epochIndx(self,epoch):
        
        # figure out where occurrences, if any, are
        indx = np.nonzero(self.epochs==epoch);
        
        # check for fail
        if len(indx) == 0:
            # just return empty
            np.array([]);
        else:
            # otherwise return the index
            return indx[0];
        
    def npvForEpoch(self,epoch):
        
        # a little bloated but easy to use epoch exists
        # but actually resolve the index twice but whatever
        if not self.epochExists(epoch):
            return np.array([]);
        else:
            return self.npvs[:,self.epochIndx(epoch)];
        
    def npvForIndx(self,indx):
        
        # get it right damnit
        if indx < 0 or indx > self.npvs.shape[1]:
            raise pyTSException('index out of bounds');
        
        # just return the npv for indx
        return self.npvs[:,np.r_[indx]];
    
    def xyzForStn(self,stns):
            
        # just snag the npv index    
        npvIndx = self.npvIndx(stns);
                
        # spit it out
        return self.npvs[npvIndx,:];
    
    def sigmaForStn(self,stns):
            
        # just snag the npv index    
        npvIndx = self.npvIndx(stns);
                
        if self.npvs_sigma == None:
            return None;
                
        # spit it out
        return self.npvs_sigma[npvIndx,:];
    
    # iterator protocol
    
    def __iter__(self):
        return self
    
    def next(self):
        
        if self.iterIndx == None:
            self.iterIndx = 0;
        
        if self.iterIndx == self.npvs.shape[1]-1:
            
            # reset iteration parameter
            self.iterIndx = None;
            
            # halt iteration
            raise StopIteration;
        else:
            self.iterIndx += 1;
            return self.npvForIndx(self.iterIndx),self.epochs[self.iterIndx];
        
    def stns(self):
        # generator for stns
        for stn in self.stn_list:
            yield stn;
            
    def spvs(self,stn):
        xyz = self.xyzForStn(stn);
        for i in range(0,xyz.shape[1]):
            if np.any(np.isnan(xyz[:,i])):
                continue;
            yield xyz[:,i],pyDate.Date(fyear=self.epochs[i]);
            
    def spvsForEpoch(self,epoch):
        
        if not self.epochExists(epoch):
            return;
        
        npv = self.npvForEpoch(epoch);
                
        for i in range(0,len(self.stn_list)):
            if np.any(np.isnan(npv[3*i:3*i+3])):
                continue;
            yield self.stn_list[i],npv[3*i:3*i+3];
            
    def spvsWithSigma(self,stn):
        xyz   = self.xyzForStn(stn);
        sigma = self.sigmaForStn(stn);
             
        for i in range(0,xyz.shape[1]):
            if np.any(np.isnan(xyz[:,i])):
                continue;
            if sigma == None:
                yield xyz[:,i],None,pyDate.Date(fyear=self.epochs[i]);
            else:
                yield xyz[:,i],sigma[:,i],pyDate.Date(fyear=self.epochs[i]);
                
    def spvsWithSigmaForEpoch(self,epoch):
        
        if not self.epochExists(epoch):
            return;
        
        npv = self.npvForEpoch(epoch);
        idx = self.epochIndx  (epoch);
                
        for i in range(0,len(self.stn_list)):
            if np.any(np.isnan(npv[3*i:3*i+3])):
                continue;
            if self.npvs_sigma == None:
                yield self.stn_list[i],npv[3*i:3*i+3,0].flatten(),None;
            else:
                yield self.stn_list[i],npv[3*i:3*i+3,0].flatten(),self.npvs_sigma[3*i:3*i+3,idx].flatten();
                
                
                
                
    