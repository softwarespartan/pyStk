

import os;
import stnInfoLib;
import pickle
import archexpl
import glob;
import re;
import sys;
import snxParse;
import HfileParser;
import PykParser;
import file_ops;
import ggParse;

class DataObjException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class DataSrcException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class DataMgrException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class DataObj():
    
    # upon initialization self.data and self.stnNameRegistry must be 
    # fully loaded and functional.  That is, XYZ data (at least) must reside
    # in self.data[stn]['xyz'] = list(x,y,z), where the station "stn" is defined
    # via forward and reverse mapping in self.stnNameRegistry.
    #
    # In general, if you just have 4char station name then initialize entry in 
    # the station name registry for igs::abcd --> abcd.
    # Notice here that IGS is the default namespace could be anything you want. 
    
    def __init__(self,expt,org,dataType,dataSrcLable,filePath):
        
        self.expt         = expt;
        self.org          = org;
        self.dataType     = dataType;
        self.dataSrcLable = dataSrcLable;
        
        self.filePath        = filePath;
        self.stnNameRegistry = None;
        
        if not os.path.isfile(self.filePath):
            raise DataObjException('file path:'+filePath+' does not exist!');  
        
        # the payload
        # data here is stored as:
        #
        # data['stn']['xyz'] = (x,y,z);
        # data['stn']['atm'] = ZND;
        # data['stn']['abc'] = ...
        # etc, etc
        #
        # this makes managing name swaps simpler to implement/coordinate
        
        self.data = dict();
                      
    def initWithStnNameRegistry(self,stnNameRegistry):
        
        self.stnNameRegistry = stnNameRegistry;
        
        return self;
        
    def initWithStnList(self,stnList):
        
        if isinstance(stnList,list):
            self.stnNameRegistry \
                = stnInfoLib.StnNameRegistry().initWithStnList(stnList);
                
        elif isinstance(stnList,dict):
            
            # make dictionary into stn list
            aList = list();
            for stnName in stnList.keys():
                aList.append(stnList[stnName]+"::"+stnName);
            
            # init the name registry
            self.stnNameRegistry \
                = stnInfoLib.StnNameRegistry().initWithStnList(aList);
        else:
            # yell at the user about it  ...
            raise DataObjException('input must be of data dataType list() or dict()!');
        
        return self;
        
    def initWithPkl(self,pklPath):
        
        if not os.path.isfile(pklPath):
            raise DataObjException('pkl path '+pklPath+' does not exist');
        
        # flag to rezip at end
        wasZipped     = False
        wasCompressed = False
        
        # check for gzip
        if pklPath[-2:] == "gz":
            file_ops.gunzip(pklPath)
            pklPath = pklPath[0:-3]
            wasZipped = True
        
        # check for unix compression
        elif pklPath[-1:] == "Z":
            file_ops.uncompress(pklPath)
            pklPath = pklPath[0:-2]
            wasCompressed = True
        
        # should put try/catch here!
        pkl_file = open(pklPath, 'rb');
        dataObj = pickle.load(pkl_file);
        pkl_file.close();
        
        if isinstance(dataObj,dict):
            if dataObj.has_key('stnLIST'):
                self.initWithStnList(dataObj['stnLIST']);
            elif dataObj.has_key('stnNameRegistry'):
                self.initWithStnNameRegistry(dataObj['stnNameRegistry']);
            else:
                raise DataObjException('could not initialize with pkl file '+pklPath);
                
        else:
            raise DataObjException('pkl data dataType not recognized!');
        
        # re-gzip the file is was .gz
        if wasZipped:
            file_ops.gzip(pklPath)
            
        # recompress the file if was .Z    
        if wasCompressed:
            file_ops.compress(pklPath)
        
        return self;
        
    def size(self):
        return self.stnNameRegistry.size()
    
    def contains(self,stnId):
        #and self.data.has_key(self.resolve(stnId))
        return self.stnNameRegistry.contains(stnId);
    
    def name(self):
        if self.dataSrcLable == "":
            return self.expt+"."+self.org+"."+self.dataType;
        else:
            return self.expt+"."+self.org+"."+self.dataType+"."+self.dataSrcLable;
    
    def resolve(self,key):
        return self.stnNameRegistry.resolve(key);
    
    def __swapStnNames(self,oldName,newName):
        
        if len(oldName) != 4 or len(newName) != 4:
            raise DataObjException('station names ' +oldName, +' and '+ newName + ' must be 4 chars');
        
        if not self.data.has_key(oldName.upper()):
            os.sys.stderr.write('stkData WARNING: '\
                                +'station name '+oldName\
                                +' does not exist in data object '\
                                +self.name()+"\n");
            return;
            #raise DataObjException('station name '+oldName+' does not exist in data object '+self.name());
            
        # first make a copy of the data under new name
        self.data[newName.upper()] = self.data[oldName.upper()];
        
        # new delete the old data
        del self.data[oldName.upper()]
    
    def updateNameMapping(self,forwardKey,stnName):
        
        if self.stnNameRegistry == None:
            raise DataObjException('stnNameRegistry is not initialize')
        
        # need to figure out what the forwardKey points 
        # at before we update it.  This is the name to replace
        # in the dataObj so don't loose it
        oldStnName = self.stnNameRegistry.resolve(forwardKey);
        
        # update mapping in the station name registry 
        try:
            self.stnNameRegistry.updateNameMapping(forwardKey,stnName);
        except Exception, e:
            print e
            #sys.stderr.write(e.value + '\n');
            return;
        
        # ok, now need to propagate the name change to the dataObj
        self.__swapStnNames(oldStnName, stnName)
        
    def __getData(self,dataKey):
        coords = dict();
        for stn in self.data.keys():
            if self.data[stn].has_key(dataKey):
                coords[stn] = self.data[stn][dataKey];
        return coords;
                
    def coordsAsDict(self):
        return self.__getData('xyz');
    
    def atmAsDict(self):
        return self.__getData('atm');
        
class DataSrc():
    
    def __init__(self,year,doy):
        
        # for debug purposes
        # to make sure that updated stn registries recognize stnIds etc
        self.newKeys = list();
        
        # iter protocol shit
        self.iterIndx = 0;
        
        # good to keep these around for name formating
        # and to make sure when making HsrcCollection that
        # all the dates match ...
        self.year = archexpl.get_norm_year_str(year);
        self.doy  = archexpl.get_norm_doy_str(doy);
        
        # list to hold dataObj objects
        self.dataObjList = list();
        
        # force no plk file
        self.NO_PKL = False;
        
    def getFileList(self, srcRoot):
        return list();
        
    def getLable(self,filePath):
        return "DataSrcLable";
        
    def addProjectSrc(self, expt, org, dataType, srcRoot):
        #
        #  Your job here is to fill up self.dataObjList
        #  with a bunch of your dataObjects from srcRoot
        #  whoes files are found by self.getFilesList(srcRoot);
        #
        # Typically, initialize data objects as
        #
        #      DataObj(expt,dataType,dataSrcLable,filePath)
        #
        # and 
        #
        #      DataObj(expt,dataType,dataSrcLable,filePath).initWithPkl(pklPath);
        #
        #
        # Example:
        #
        #    for file in self.getFilesList():
        #        dataSrcLable = self.getLable(filePath);
        #        self.dataObjList                                                      \
        #            .append(                                                          \
        #                    DataObj( expt, org, dataType, dataSrcLable, filePath )    \
        #                        .initWithPkl( pklPath ));                             \
        #
        return self;
        
        # iterator protocol
    
    def __iter__(self):
        return self
    
    # iterator protocol
    
    def next(self):
        
        if self.iterIndx > len(self.dataObjList)-1:     
            # reset iteration parameters
            self.iterIndx = 0;
            
            # halt iteration
            raise StopIteration;
        else:
            self.iterIndx += 1;
            return self.dataObjList[self.iterIndx-1];

    def size(self):
        return len(self.dataObjList);
    
    def __len__(self):
        return len(self.dataObjList);
    
    def count(self):
        return len(self.dataObjList);
    
    def __getitem__(self,key):
        
        if not isinstance(key,int):
            try:
                # might be np.array etc
                key = int(key)
            except:
                raise KeyError;
        
        if key < 0 or key > self.size():
            raise IndexError;
        
        return self.dataObjList[key];  
    
    def getStnIdSet(self):
        
        stnIdSet = set();
        for dataObj in self:
            for stnId in dataObj.stnNameRegistry:
                if stnId.startswith('pub:'):
                    #print stnId,
                    stnId = 'igs:'+stnId[4:];
                    #print '-->',stnId
                stnIdSet.add(stnId);
        return stnIdSet
    
    def getStnNamesAsList(self):
        
        stnIdList = list();
        for dataObj in self:
            for stnId in dataObj.stnNameRegistry:
                stnIdList.append(dataObj.stnNameRegistry.resolve(stnId));
        return stnIdList

    def getForwardAliasMap(self):
        
        # The idea here is to generate a list of names that each stnId points to
        #
        #  So for example, suppose igs::algo were in three solutions
        #  In two solutions it was processed as 'algo' but in the third
        #  solution it was renamed to _xud to avoid conflict with "chi::algo".
        #
        #  Therefore the forwardAliasMap for igs::algo would look like
        #
        #    igs::algo --> set('algo','algo','_xud',) = ('algo', '_xud')
        #
        #  and likewise
        #
        #    chi::algo --> set('algo') = ('algo')
        #
        #  Notice here that the Set gives us unique entries. 
        #
        #
        #  The reverseAliasMap will provide the look up to identify stations
        #  with different names that actually map to the same station id.
        #
        #  For example, the reverse alias map for station name algo will be
        #
        #    algo --> (igs::algo, chi::algo)
        # 
        #  and likewise
        # 
        #    _xud --> (igs::algo)
        #
        #
        #  Clearly, this is a name conflict that would need to be resolved.
        
        
        # list of all uniq. stnId's 
        # note pub: --> igs:
        stnIdSet = self.getStnIdSet();
        
        # init a map s.t. idMap[stnId] = set()
        # where set contains all aliases for stnId
        idMap = dict();
        for stnId in stnIdSet:
            
            # init the alias set for this stnId
            idMap[stnId] = set();
            
            # now go through each dataObj 
            for dataObj in self:
                
                # get alias according to this registry
                alias  = dataObj.stnNameRegistry.resolve(stnId);
                
                # if the stnId resolution was successful 
                # then add it to alias set
                if alias !=None:
                    idMap[stnId].add(alias);
                
        return idMap
        
    def getReverseAliasMap(self):
        
        # get the forward map, notice that
        # if we resolve forward alias mapping conflicts
        # then might not have any revers alias mapping conflicts
        # so always good to compile forward map first
        forwardMap = self.getForwardAliasMap();
        
        #print 'forwardMap size:',len(forwardMap.keys())
        
        # init unique set of alias'
        aliasSet = set();
        for stnId in forwardMap.keys():
            aliasSet = aliasSet.union(forwardMap[stnId]);
            
        #print 'aliasSet size:',len(aliasSet);
        reverseMap = dict();
        for alias in aliasSet:
            reverseMap[alias] = set();
            
            for dataObj in self:
                stnId = dataObj.stnNameRegistry.resolve(alias);
                if stnId != None:
                    # make sure that igs:: and pub:: map the same
                    if stnId.startswith('pub:'):
                        stnId='igs:'+stnId[4:];
                    reverseMap[alias].add(stnId);
                    
        return reverseMap    

    def resolveForwardAliasConflicts(self):
        
        # OK, for forward conflicts we have uniq. stnId's
        # that point to the different keys.  
        #
        # For example:
        #
        #   igs::algo  --> algo  (proj #1)
        # -----------------------------------
        #   chi::algo  --> algo  (proj #2)
        # -----------------------------------
        #   igs::algo  --> algo  (proj #3)
        #   chi::algo  --> xjsl  
        # -----------------------------------
        #
        # So here the forward mappings:
        #
        #   chi::algo --> algo  (h-file #2)
        #   chi::algo --> xjsl  (h-file #3)
        #
        # are conflicting but are same station.
        #
        #  The resolution is to remap *ALL* chi::algo keys to 
        #  a single non conflicting key.  For example
        #
        #            chi::algo --> sqki
        #
        #  Thus the final forward map for projects 1 - 3 for stnId's
        #  igs::algo and chi::algo is:
        #
        #   igs::algo  --> algo  (proj #1)
        # -----------------------------------
        #   chi::algo  --> sq_i  (proj #2)
        # -----------------------------------
        #   igs::algo  --> algo  (proj #3)
        #   chi::algo  --> sq_i  
        # ----------------------------------- 
        #
        #  In words, all occurances of ALGO in h-file #2 will be
        #  be replaced with SQ_I.  And, all accurances of XJSL in
        #  h-file #3 will be replaced with SQ_I.
        
        # get the forward alias mapping
        fMap = self.getForwardAliasMap();
        
        # now loop throgh and check reference count
        for id in fMap.keys():
            
            # if the reference count > 1 = conflict!
            if len(fMap[id])>1:
                
                # print out some debug shit
                print 'FORWARD CONFLICT: ',id,'-->',
                for alias in fMap[id]:
                    print alias+',',
                print
                
                # make up new name
                newName = stnInfoLib.StnNameRegistry().getRandom4Chars();
                print "RESOLUTION:", id,'-->',newName;
                
                # assign this to every h file in the collection
                for dataObj in self:
                    if dataObj.contains(id):
                        print 'RESOLUTION: updating  '+dataObj.name();
                        dataObj.updateNameMapping(id,newName);
                        
    def resolveReverseAliasConflicts(self):
        
        # Reverse alias mapping conflicts are about the same difficulty
        #
        # An occurance of a reverse alias conflict is such that two 
        # indepenant projects process a station of the same name but resides
        # in two different name spaces.  For example:
        #
        #   igs::kouc  --> kouc   (proj #1)
        # -----------------------------------
        #   swp::kouc  --> kouc   (proj #2)
        # -----------------------------------
        #   bra::kouc  --> kouc   (proj #3)
        #
        # So here the key 'kouc' has reverse mapping:
        #
        #   kouc --> igs::kouc, swp::kouc, bra::kouc
        #
        # NOTE: The example in the resolveForwardAliasConflicts function above
        #       had a reverseAliasConflict but was resolved during forward-resolution!
        #       It's a bit subtle but make sure you understand!!!  (or am i confused?)
        #       Anyway ...
        #
        # The resolution here is remap N-1 stnId's to new keys.  Notice that if 
        # we reassign swp::kouc and bra::kouc that the mapping igs::kouc --> kouc
        # is no longer in conflict.  Thus 3 reverse mapping conflicts require 2 updates. 
        #
        # Thus the final mapping for kouc after reverseAliasConflict resolution is:
        #
        #   igs::kouc  --> kouc   (proj #1)
        # -----------------------------------
        #   swp::kouc  --> pw_k   (proj #2)
        # -----------------------------------
        #   bra::kouc  --> mc_v   (proj #3)
        #
        #
        
        # compute the reverse alias mapping
        rMap = self.getReverseAliasMap();

        for  alias in rMap.keys():

            # if the reference count for this alias is
            # greater than 1 we have a conflict
            if len(rMap[alias])>1:
                
                # print out debug shit for user
                print 'REVERSE CONFLICT:  ',alias,'-->',
                for stnId in rMap[alias]:
                    print stnId+',',
                print
                
                for stnId in rMap[alias]:
                    
                    newName = stnInfoLib.StnNameRegistry().getRandom4Chars();
                    self.newKeys.append(newName)
                    print 'RESOLUTION:',stnId,'-->',newName
                    
                    for dataObj in self:
                        if dataObj.contains(stnId):
                            print 'RESOLUTION: updating  '+dataObj.name();
                            dataObj.updateNameMapping(stnId,newName);

    def resolveNameConflicts(self):
        self.resolveForwardAliasConflicts();
        self.resolveReverseAliasConflicts();
        
    def getReferenceIndx(self,key):
        
        i=0;
        indx = list();
        for dataObj in self:
            if dataObj.contains(key.lower()):
                indx.append(i);
            i = i+1;
        return indx;
            
    def getReferenceCount(self,key):
        
        return len(self.getReferenceIndx(key));
          
    def getReferences(self,key):
        l = list();
        for indx in self.getReferenceIndx(key):
                l.append(self.dataObjList[indx]);
        return l;                

class DataMgr(DataSrc):
    
    # a matter of perpective really ...
    
    def __init__(self, year, doy):
        DataSrc.__init__(self, year, doy);
        
    def addData(self,dataSrcObj):
        
        # make sure that the dates are the same
        if dataSrcObj.year != self.year \
            or dataSrcObj.doy != self.doy :
            raise DataMgrException('Can not add data object with different dates');

        # just copy over each data object from the source
        for dataObj in dataSrcObj:
            self.dataObjList.append(dataObj);
            
class SnxDataObj(DataObj):
    
    def __init__(self, expt, org, dataType, dataSrcLable, filePath):
        
        # call super class constructor for  data object
        DataObj.__init__( self, expt, org, dataType, dataSrcLable, filePath );
        
        # init sinex parser for current sinex file   
        snxParser = snxParse.snxFileParser(filePath).parse();
        
        # organize the coordinates for the data object
        for stn in snxParser:
            
            # initialize the station if we haven't seen it before
            if not self.data.has_key(stn):
                self.data[stn] = dict();

            # get the snx data from the parser
            stnData = snxParser.get(stn);
            
            # pack the data away
            self.data[stn]['xyz'] = [stnData.X,stnData.Y,stnData.Z];
            
class SnxSrc(DataSrc):
    
    def __init__(self, year, doy):
    
        DataSrc.__init__(self, year, doy);
        
    def getFileList(self, srcRoot, org = ""):
        
        # pattern that defines a file
        fileDef = org+"*.snx*"
        
        # generate date specific path from srcRoot + org
        #  option 1:  year/doy ...
        #  option 2:  year/doy/org ...
        hpathYearDoy = os.path.join(os.path.join(srcRoot,self.year),self.doy);
        hpathYearDoyOrg = os.path.join(hpathYearDoy,org);
        
        # choose which path to use
        if org == "" or org == None:
            
            # solutions are in bare year/doy path
            hpath = hpathYearDoy;
        
        elif os.path.isdir(hpathYearDoyOrg):
            
            # ok, now the solutions will be in year/doy/org path
            hpath = hpathYearDoyOrg;
        
        else:
            # dflt
            hpath = hpathYearDoy;
        
        # if the path does not exist then we're done
        if not os.path.isdir(hpath):
            return list();
            
        # construct a path for each network
        net_path = os.path.join(hpath,'n[0-9]*');
        
        # look for .snx file in each network direcory
        path = os.path.join(net_path,fileDef)
        
        # look for the files in year/doy/netid/
        project_list = glob.glob(path);
        
        # if we didn't find any
        if len(project_list) == 0:
            
            # then just look in the original hpath
            path= os.path.join(hpath,fileDef);
            
            # look for files in year/doy dir     
            project_list = glob.glob(path);
         
        # and we're done ...   
        return project_list;
    
    def getLable(self,filePath):
        regex = re.compile("(n\d+)")
        match = regex.findall(os.path.dirname(filePath))
        if len(match) == 0:
            return "";
        if len(match) == 1:
            return match[0]
    
    def getPkl(self,filePath,org=None):
        
        # descript the file we're looking for
        fileDef = "*.pkl*";
        
        # add the org if provided 
        if org != None:
            fileDef = org+fileDef;
        
        # create the directory path which should have pkl file
        fileDir = os.path.split(filePath)[0];
            
        # next check for pkl file for this h file
        pklFileList = glob.glob(os.path.join(fileDir,fileDef));
        
        # check that hFile has exactly 1 pkl file
        if len(pklFileList) == 1:
            return pklFileList[0]
        else:
            return None;
    
    def addProjectSrc(self, expt, org, dataType, srcRoot):
        
        # for each sinex file in the solution root ...
        for snxFile in self.getFileList(srcRoot, org):
            
            # get lable associated with snxFile path
            dataSrcLable = self.getLable(snxFile);
            
            # get the pkl file associated with this sinex file
            pklPath = self.getPkl(snxFile);
            
            # create a sinex object 
            snxObj = SnxDataObj(expt, org, dataType, dataSrcLable, snxFile);
            
            # if there is no pkl file then move along
            if pklPath == None :
                #sys.stderr.write('snxRoot: '+snxFile+' does not have appropriate pkl associated with it!');
                #continue;
                stnList = list();
                for s in snxObj.data.keys():
                    stnList.append('igs::'+s);
                    
                self.dataObjList.append(snxObj.initWithStnList(stnList));
            else:
                # load the sinex file
                self.dataObjList.append(SnxDataObj(expt, org, dataType, dataSrcLable, snxFile).initWithPkl(pklPath));
                        
        return self;
    
class HfileDataObj(DataObj):
    
    def __init__(self, expt, org, dataType, dataSrcLable, filePath):
        
        # call super class constructor for  data object
        DataObj.__init__( self, expt, org, dataType, dataSrcLable, filePath );
        
        # init sinex parser for current  hfile   
        hfileParser = HfileParser.HfileParser(filePath).parse(dataType);
        
        # organize the coordinates for the data object
        for stn in hfileParser:
            
            # initialize the station if we haven't seen it before
            if not self.data.has_key(stn):
                self.data[stn] = dict();

            # get the snx data from the parser
            stnData = hfileParser.get(stn);
            
            # pack the data away
            self.data[stn]['xyz'] = [stnData.X,stnData.Y,stnData.Z];
                    
class HfileSrc(DataSrc):
    
    def __init__(self, year, doy):
        DataSrc.__init__(self, year, doy)
        
    def getFileList(self, srcRoot, org = ""):
        
        # pattern that defines a file
        fileDef = 'h*.[0-9][0-9][0-9][0-9][0-9]*'
        
        # generate date specific path from srcRoot + org
        #  option 1:  year/doy ...
        #  option 2:  year/doy/org ...
        hpathYearDoy    = os.path.join(os.path.join(srcRoot,self.year),self.doy);
        hpathYearDoyOrg = os.path.join(hpathYearDoy,org);
                
        # choose which path to use
        if org == "" or org == None:
            
            # solutions are in bare year/doy path
            hpath = hpathYearDoy;
                    
        elif os.path.isdir(hpathYearDoyOrg):
            
            # ok, now the solutions will be in year/doy/org path
            hpath = hpathYearDoyOrg;
                     
        else:
            # dflt
            hpath = hpathYearDoy;
        
        # if the path does not exist then we're done
        if not os.path.isdir(hpath):
            return list();
            
        # construct a path for each network
        net_path = os.path.join(hpath,'n[0-9]*');
        
        # look for .snx file in each network direcory
        path = os.path.join(net_path,fileDef)
        
        # look for the files in year/doy/netid/
        project_list = glob.glob(path);
        
        # if we didn't find any
        if len(project_list) == 0:
            
            # then just look in the original hpath
            path= os.path.join(hpath,fileDef);
            
            # look for files in year/doy dir     
            project_list = glob.glob(path);
                     
        # and we're done ...   
        return project_list
       
    def getLable(self,filePath):
        regex = re.compile("(n\d+)")
        match = regex.findall(os.path.dirname(filePath))
        if len(match) == 0:
            return "";
        if len(match) == 1:
            return match[0]        
    
    def getPkl(self,filePath):
        
        fileDir = os.path.split(filePath)[0];
            
        # next check for pkl file for this h file
        pklFileList = glob.glob(os.path.join(fileDir,'*.pkl'));
        
        # check that hFile has exactly 1 pkl file
        if len(pklFileList) == 1:
            return pklFileList[0]
        else:
            sys.stderr.write('hfileRoot: '+filePath+' does not have appropriate pkl associated with it!');
            return None;
            
    def addProjectSrc(self, expt, org, dataType, srcRoot):
        
        # for each sinex file in the solution root ...
        for hfile in self.getFileList(srcRoot, org):
            
            # get the pkl file associated with this sinex file
            pklPath = self.getPkl(hfile);
            
            # if there is no pkl file then move along
            if pklPath == None:
                continue;
            
            # get lable associated with snxFile path
            dataSrcLable = self.getLable(hfile);
            
            print "adding file: ",hfile
            
            # load the sinex file
            self.dataObjList.append(HfileDataObj(expt, org, dataType, dataSrcLable, hfile)
                                    .initWithPkl(pklPath));
                        
        return self;
    
class PykDataObj(DataObj):
    
    def __init__(self, expt, org, dataType, dataSrcLable, filePath):
        
        # call super class constructor for  data object
        DataObj.__init__( self, expt, org, dataType, dataSrcLable, filePath );
        
        # init sinex parser for current  hfile   
        parser = PykParser.PykParser(filePath).parse(dataType);
        
        # need to initialize a station name registry for these folks
        self.stnNameRegistry \
            = stnInfoLib.StnNameRegistry().initWithStnList(parser.stationDict.keys());
        
        # organize the coordinates for the data object
        for stnId in parser:
            
            # igs::zimm -> zimm
            stnName = self.stnNameRegistry.resolve(stnId)
            
            # initialize the station if we haven't seen it before
            if not self.data.has_key(stnName):
                self.data[stnName] = dict();

            # get the snx data from the parser
            stnData = parser.get(stnId);
            
            # pack the data away
            self.data[stnName]['xyz'] = [stnData.X,stnData.Y,stnData.Z];
            
class PykSrc(DataSrc):
    
    def __init__(self, year, doy):
        DataSrc.__init__(self, year, doy)
        
    def getFileList(self, srcRoot, org = ""):
        
        # pattern that defines an hfile
        fileDef = org+self.year+self.doy+'.pyk*';
                
        # generate date specific path from srcRoot
        hpath = srcRoot;
    
        # if the path does not exist then we're done
        if not os.path.isdir(hpath):
            return list();
            
        # construct a path for each network
        net_path = os.path.join(hpath,'n[0-9]*');
        
        # look for file in each network direcory
        path = os.path.join(net_path,fileDef)
        
        # look for the files in year/doy/netid/
        project_list = glob.glob(path);
        
        # if we didn't find any
        if len(project_list) == 0:
            
            # then just look in the original hpath
            path= os.path.join(hpath,fileDef);
            
            # look for files in year/doy dir     
            project_list = glob.glob(path);
            
        return project_list;
       
    def getLable(self,filePath):
        
        regex = re.compile("(n\d+)")
        
        match = regex.findall(os.path.dirname(filePath))
        
        if len(match) == 0:
            return "";
        if len(match) == 1:
            return match[0]
            
    def addProjectSrc(self, expt, org, dataType, srcRoot):
        
        # for each sinex file in the solution root ...
        for file in self.getFileList(srcRoot, org):
            
            # get lable associated with snxFile path
            dataSrcLable = self.getLable(file);
            
            # load the sinex file
            self.dataObjList.append(PykDataObj(expt, org, dataType, dataSrcLable, file));
                        
        return self;    
    
class PrtDataObj(DataObj):
    
    def __init__(self, expt, org, dataType, dataSrcLable, filePath):
        
        # call super class constructor for  data object
        DataObj.__init__( self, expt, org, dataType, dataSrcLable, filePath );
        
        # init sinex parser for current  hfile   
        parser = ggParse.prtFileParser(filePath);
        
        # parse the prt file
        parser.parse();
        
        # set the xyz coordinates for the stations
        self.data = parser.stationCoordinates;
            
class PrtSrc(DataSrc):
    
    def __init__(self, year, doy):
        DataSrc.__init__(self, year, doy)
        
    def getFileList(self, srcRoot, org = ""):
        
        # pattern that defines an hfile
        fileDef = 'gkglbk'+self.doy+self.year[-2:]+'.prt*';
        
        # generate date specific path from srcRoot + org
        #  option 1:  year/doy ...
        #  option 2:  year/doy/org ...
        hpathYearDoy    = os.path.join(os.path.join(srcRoot,self.year),self.doy);
        hpathYearDoyOrg = os.path.join(hpathYearDoy,org);
        
        # choose which path to use
        if org == "" or org == None:
            
            # solutions are in bare year/doy path
            hpath = hpathYearDoy;
                    
        elif os.path.isdir(hpathYearDoyOrg):
            
            # ok, now the solutions will be in year/doy/org path
            hpath = hpathYearDoyOrg;
                     
        else:
            # dflt
            hpath = hpathYearDoy;
        
        # if the path does not exist then we're done
        if not os.path.isdir(hpath):
            return list();
            
        # construct a path for each network
        net_path = os.path.join(hpath,'n[0-9]*');
        
        # look for .snx file in each network direcory
        path = os.path.join(net_path,fileDef)
        
        # look for the files in year/doy/netid/
        project_list = glob.glob(path);
        
        # if we didn't find any
        if len(project_list) == 0:
            
            # then just look in the original hpath
            path= os.path.join(hpath,fileDef);
            
            # look for files in year/doy dir     
            project_list = glob.glob(path);
                     
        # and we're done ...   
        return project_list
       
    def getLable(self,filePath):
        
        regex = re.compile("(n\d+)")
        
        match = regex.findall(os.path.dirname(filePath))
        
        if len(match) == 0:
            return "";
        if len(match) == 1:
            return match[0]
            
    def getPkl(self,filePath,org=None):
        
        # descript the file we're looking for
        fileDef = "*.pkl*";
        
        # add the org if provided 
        if org != None:
            fileDef = org+fileDef;
        
        # create the directory path which should have pkl file
        fileDir = os.path.split(filePath)[0];
            
        # next check for pkl file for this h file
        pklFileList = glob.glob(os.path.join(fileDir,fileDef));
        
        # check that hFile has exactly 1 pkl file
        if len(pklFileList) == 1:
            return pklFileList[0]
        else:
            return None;
            
    def addProjectSrc(self, expt, org, dataType, srcRoot):
        
        # for each sinex file in the solution root ...
        for file in self.getFileList(srcRoot, org):
                        
            # get the pkl file associated with this sinex file
            pklPath = self.getPkl(file); 
                        
            # if there is no pkl file then move along
            if pklPath == None:
                os.sys.stderr.write('no pkl for file '+file+'\n');
                continue;            
                        
            # get lable associated with snxFile path
            dataSrcLable = self.getLable(file);
            
            print "adding file: ",file
            
            # load the sinex file
            self.dataObjList                                                    \
                .append(PrtDataObj(expt, org, dataType, dataSrcLable, file)     \
                        .initWithPkl(pklPath)                                   \
            );
                        
        return self;    