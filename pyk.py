#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import os
import sys, traceback
import getopt
import archexpl;
import globk_ops;
import hfileLib;
import pykHelper;
import file_ops;
import traceback;
import stkData;

class PykException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

help_short_arg        = "h"
help_long_arg         = "help="

outdir_short_arg      = "o:"
outdir_long_arg       = "outdir="

file_short_arg        = "f:"
file_long_arg         = "file="

compress_short_arg   = "c"
compress_long_arg    = "compress"

year_short_arg       = "y:"
year_long_arg        = "year="

doy_short_arg        = "d:"
doy_long_arg         = "doy="

frame_long_arg       = "rf="

# shlep'em all together
short_args = help_short_arg    + \
             outdir_short_arg  + \
             file_short_arg    + \
             year_short_arg    + \
             doy_short_arg     + \
             compress_short_arg;
             
long_args = []
long_args.append(help_long_arg)
long_args.append(outdir_long_arg)
long_args.append(file_long_arg);
long_args.append(compress_long_arg);
long_args.append(year_long_arg);
long_args.append(doy_long_arg);
long_args.append(frame_long_arg);

def usage():
    
    print
    print "USAGE:  pyk -y yyyy -d ddd -f /path/to/pyk.config/file"
    print 
    print "options: "
    print
    print "-f, --file=         path to config file"
    print "-e, --expt=         prefix for solution file name"
    print "-o, --outdir=       specify output directory for .pyk file.  DEFAULT: pwd"
    print "-c, --compress      UNIX compress output file "
    print "  , --rf=           Impose reference frame on solution (itrf08/itrf05)"
    print "-h, --help          show this help message and exit"
    print 

def getInputArgs():
    
    # defaults
    file            = None;
    expt            = '';
    outdir          = None;
    shouldCompress  = False;
    
    year = None;
    doy  = None;
    
    refFrame = None;
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], short_args, long_args)
        
    except getopt.GetoptError, err:
        
        # print help information and exit:
        print str(err)
        usage()
        sys.exit(2)
        
    if len(opts)==0:
        usage();
        sys.exit();
        
    for option,arg in opts:
        
        # help arg actions
        if option in ("-h","--help"):
                usage()
                sys.exit(1)
                
        if option in ("-y","--year"):
            
            try:
                year = int(arg)
            except:
                raise PykException('Could not parse year argument')
                
            if year < 1980:
                raise PykException('Invalid year');
                
        elif option in ("-d","--doy"):
            try:
                doy = int(arg)
            except:
                raise PykException('Could not parse doy argument')
                
            if doy < 1 or doy > 366:
                raise PykException('Invalid day-of-year');
                
        elif option in ("-e","--expt"):
            expt = arg;
        
        elif option in ("-o","--outdir"):
            
            # assign outdir
            outdir = arg
            
            # cross platform junk
            outdir = os.path.normcase(outdir)
            outdir = os.path.normpath(outdir)
            
            # verify 
            if not os.path.isdir(outdir):
                raise PykException("outdir "+outdir+" does not exist")
                
        elif option in ("-f","--file"):
            file = arg;
            if not os.path.isfile(file):
                raise PykException("file "+file+" does not exist");
            
        elif option in ("-c","--compress"):
            shouldCompress = True;
       
        elif option in ("--rf"):
            if arg.lower() in ('itrf08','itrf2008','itrf05' 'itrf2005'):
                refFrame = arg.lower();
                
            
    return (year,doy,file,expt,outdir,shouldCompress,refFrame)

def ckInputArgs(year,doy,file,expt,outdir,shouldCompress,refFrame):
    
    # gotta have a date to know where to get data
    if year == None or doy == None:
        raise PykException('Must specify year and doy using -y and -d options');
    
    # no config = deal breaker!
    if file == None:
        raise PykException('Must specify path to pyk.config file using -f option');
        
    return year,doy,file,expt,outdir,shouldCompress,refFrame

def getOutFileName(year,doy,expt=''):
    return expt+str(year)+archexpl.get_norm_doy_str(doy)+'.pyk';

def initPykConfig(file):
    
    return globk_ops.read_globk_config(file);
    
def initDataSrc(year,doy,pykConfig):
    
    # init src var
    src = None;
    
    # Mgr of data sources
    data = stkData.DataMgr(year,doy);
                
    # for each DataSrc listed in .config file ...    
    for solnSrc in pykConfig['DataSrc']:
    
        # parse the line from .config file
        (dataSrcTag, srcType, expt, org, dataType, srcRoot) = solnSrc.split('::');
                
        # initialize the DataSrc object according to srcType
        if srcType.lower() == 'snxsrc':
            src = stkData.SnxSrc(year,doy);
        elif srcType.lower() == 'hfilesrc':
            src = stkData.HfileSrc(year,doy);
        elif srcType.lower() == 'pyksrc':
            src = stkData.PykSrc(year,doy);
        else:
            raise PykException(srcType + " is unrecognized data source type for Pyk");
                    
        try:
            # add any data here to collection
            data.addData(src.addProjectSrc(expt, org, dataType, srcRoot));
        except Exception, e:
            print e;
            traceback.print_exc(file=sys.stdout)
            sys.stderr.write('Error importing data files from data source type: '+srcType.lower()+": "+srcRoot+'\n');
        
    # finally, fix any name collisions
    if pykConfig['shouldResolveNameConflicts']:
        data.resolveNameConflicts();
    
    return data;
    
if __name__ == "__main__":
    
    # round'em up
    year,doy,file,expt,outdir,shouldCompress,refFrame = getInputArgs();
    
    # check'em out
    year,doy,file,expt,outdir,shouldCompress,refFrame \
        = ckInputArgs(year,doy,file,expt,outdir,shouldCompress,refFrame);

    # read the configuration file
    pykConfig = initPykConfig(file);
    
    # initialze all the data in each DataSrc
    data = initDataSrc(year,doy,pykConfig);
    
    # make sure some data was initialized 
    if data.size() == 0:
        sys.stderr.write('no data found for '+str(year)+" "+str(doy)+'\n');
        sys.exit(1);
    
    # get those data ready
    helper = pykHelper.helper(data); 
    
    try:
        # do the deed
        helper.stk();
    except Exception as e:
        traceback.print_exc()
        sys.exit(2)
        
    # make sure nothing sneaky happend (?)
    if not helper.isValid():
        sys.exit(2);
        
    # figure out file prefix
    if expt == '' and pykConfig.has_key('expt'):
        expt = pykConfig['expt'];
        
    # generate a name for output
    outFileName = getOutFileName(year,doy,expt)
    
    if outdir != None:
        outFilePath = os.path.join(outdir,outFileName);
    elif pykConfig['solutionROOT'] != None:
        outFilePath = os.path.join(pykConfig['solutionROOT'],outFileName);
    else:
        outFilePath = os.path.join(os.path.normpath("."),outFileName);
         
        
    # initialize the output file
    fid = open(outFilePath,'w');
    
    # generate the output 
    helper.printSummary(fid);
    
    if refFrame != None:
        helper.imposeRF(refFrame, year, doy, fid);
    
    helper.printCoords(fid);
    
    # close the output files, we're done here
    fid.close();
    
    # lastly, compress the output file if desired
    if shouldCompress:
        file_ops.compress(outFilePath);
