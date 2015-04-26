

import numpy as np;
import scipy.stats as stats;

def ipr(dv,percentile=None):
    
    # compute the (I)nter(P)ercentile(R)ange
    
    if percentile == None:
        percentile = 50.0;
    
    u_ip = stats.scoreatpercentile(dv[~np.isnan(dv)], 100 - percentile/2);
    l_ip = stats.scoreatpercentile(dv[~np.isnan(dv)],       percentile/2);    
    
    return  u_ip - l_ip;

class weights():
    
    def __init__(self):
        
        self.percentile = 50;
        self.limit      = 1.65;
        self.reweight   = 1.5;
        
        # if over by greater than maxLimit sigmas 
        # then weight = 0 (stop loss) 
        #
        #  That is, 
        # 
        #     if residual_i > ipRange*maxLimit, then w_i = 0;
        #
        self.maxLimit   = 10;
        
        # a priori aprioriSigma
        self.aprioriSigma = None;
            

    def computeWeights(self,dv):

        weights = np.ones_like(dv);
        
        # compute the inter-percentile range
        if self.aprioriSigma == None:
            effectiveSigma = ipr(dv,self.percentile);
                
            #print "PYWEIGHTS: sigma_ipr = ",effectiveSigma;
                
        else:
            effectiveSigma = self.aprioriSigma;
            #print "PYWEIGHTS: sigma_apriori = ",effectiveSigma;      
            
        if effectiveSigma < 1e-4:
            return weights;  
                        
        # figure out where data over the limit is
        indx = np.abs(dv)>effectiveSigma*self.limit;
        
        # figure out where the outliers are 
        indxMax = np.abs(dv)> effectiveSigma * self.maxLimit;
        
        # no measurment is overlimit
        if np.nonzero(indx)[0].size == 0:
            return weights;
        
        # compute the weight factor for each over limit dv
        weightFactor = ( ( np.abs( dv[indx] ) / effectiveSigma) ) - self.limit;
    
        #compute the weights
        weights[indx] =  self.reweight ** -weightFactor;
        
        # set outlier weights to zero
        if np.nonzero(indxMax)[0].size > 0: 
            weights[indxMax] = 0; 
        
        return weights