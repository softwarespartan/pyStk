
import referenceFrames;
import pyStk;
import re;
import numpy as np;

class pyRFException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class pyRF():
    
    def __init__(self,refFrame):
        
        self.refStnList = list();
        self.refData = np.array([]);
        
        # get the reference frame definition
        if refFrame.lower() in ('itrf08','itrf2008'):
            self.refFrame = referenceFrames.itrf2008;
        elif refFrame.lower() in ('itrf05','itrf2005'):
            self.refFrame = referenceFrames.itrf2005;
        else:
            raise pyRFException('Reference frame '+refFrame+' not recognized');
        
        # save the id 
        self.refFrameId = refFrame;   
        
        # to be determined ...
        self.stnList = [];
        self.npv = np.array([]);
        self.nvv = np.array([]);
        self.refEpoch = np.array([]);
        
        self.__parse();
        
    def __parse(self):
                
        # parse the line on spaces
        lineParts = re.split('\s+',self.refFrame);
        del lineParts[0];

        for i in range(0,len(lineParts)/8):
            
            # cast the line parts
            stnId = lineParts[i*8];
            
            x = float(lineParts[i*8+1]); 
            y = float(lineParts[i*8+2]);
            z = float(lineParts[i*8+3]);
            
            vx = float(lineParts[i*8+4]);
            vy = float(lineParts[i*8+5]);
            vz = float(lineParts[i*8+6]);
            
            refEpoch = float(lineParts[i*8+7]);
        
            # add the stnid to reference list
            self.refStnList.append(stnId.lower());
            
            # add the data
            if self.refData.size == 0:
                self.refData = np.array([x,y,z,vx,vy,vz,refEpoch]);
            else:
                self.refData = np.vstack((self.refData,np.array([x,y,z,vx,vy,vz,refEpoch])));
                
    def initForStnList(self,stnList):
        
        self.stnList = stnList;
        
        # mem alloc
        numStns       = len(stnList);
        self.npv      = np.zeros(3*numStns,dtype='float')*np.nan;
        self.nvv      = np.zeros(3*numStns,dtype='float')*np.nan;
        self.refEpoch = np.zeros(3*numStns,dtype='float')*np.nan;
        
        for i in range(0,len(stnList)):
            
            stn = stnList[i];
            
            if stn in self.refStnList:
                
                refIndx = self.refStnList.index(stn);
                indx    = 3*i;
                
                self.npv[indx:indx+3]      = self.refData[refIndx,[0,1,2]];
                self.nvv[indx:indx+3]      = self.refData[refIndx,[3,4,5]];
                self.refEpoch[indx:indx+3] = self.refData[refIndx,[6]].repeat(3);
                                
        return self;
    
    def init(self):
    
        self.initForStnList(self.refStnList);
        return self;
    
    def npvForEpoch(self,fyear):
        
        dt  = fyear - self.refEpoch;
        npv = self.npv + (self.nvv * dt);
        npv = npv.transpose().reshape(npv.size,1);
        
        return npv;
    
    def transformForEpoch(self,npv,fyear):
        
        npvTarget    = self.npvForEpoch(fyear);
        T,npvT,stats = pyStk.helmert(npv, npvTarget);
        return T,npvT,stats
    
    def transformWithScaleForEpoch(self,npv,fyear):
        
        npvTarget    = self.npvForEpoch(fyear);
        T,npvT,stats = pyStk.helmert(npv, npvTarget,withScale='yes');
        return T,npvT,stats
            