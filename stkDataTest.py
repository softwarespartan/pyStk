
import stkData;

src     = "../data/osuGLOBAL/napeos/densification/solutions/";
snxfile = "../data/osuGLOBAL/napeos/orbit/solutions/1998/050/os609454.snx.Z";
pklPath = "../data/osuGLOBAL/napeos/orbit/solutions/1998/050/glbl1998050.pkl"

#snxData = stkData.SnxDataObj('glbl', 'snx', 'glbl', snxfile).initWithPkl(pklPath);
#coords = snxData.coordsAsDict();
#
#for stn in coords.keys():
#    print snxData.name(), stn, coords[stn];

snxSrc = stkData.SnxSrc(1998,50).addProjectSrc('glbl', 'snx', src);

for obj in snxSrc:
    print obj.name();
    
#snxSrc.resolveNameConflicts();
#print len(snxSrc.getStnIdSet());

