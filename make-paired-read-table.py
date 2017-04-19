import glob
import sys
from optparse import OptionParser
###############################################################################
USAGE = """
make-paired-read-table.py  --in <dir>   

will make read-pair table of paired reads.  Then you should sort the table and convert
to fq.gz

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inDir', help = 'input dir to process')


(options, args) = parser.parse_args()

if options.inDir is None:
    parser.error('inDir not given')

###############################################################################

if options.inDir[-1] != '/':
    options.inDir += '/'


isPaired = {}
pairedListFile = options.inDir +'info.parsed.paired'
inFile = open(pairedListFile,'r')
for line in inFile:
    line = line.rstrip()
    isPaired[line] = 1
inFile.close()
print 'Read in list of %i that are paired' % len(isPaired)


tiToName = {}
infoFile = options.inDir +'info.parsed'
inFile = open(infoFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    ti = line[0]
    cName = line[1]
    cDir = line[2]
    if cName in isPaired:
        if cDir == 'REVERSE':
            n = cName + ' 2'
            r = 1
        elif cDir == 'FORWARD':
            n = cName + ' 1'
            r = 0
        else:
            print 'What dir?'
            print line
            sys.exit()
        tiToName[ti] = (n,r)
inFile.close()
print 'Read in  %i that are ready to go' % len(tiToName)

read1TableName = options.inDir +'seq_qual.R1'
read2TableName = options.inDir +'seq_qual.R2'

outFiles = [0,0]
outFiles[0] = open(read1TableName,'w')
outFiles[1] = open(read2TableName,'w')
seqQualFileName = options.inDir +'seq_qual.parsed'
inFile = open(seqQualFileName,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    ti = line[0]
    seq = line[1]
    qual = line[2]
    if ti in tiToName:
        n = tiToName[ti][0]
        fn = tiToName[ti][1]
        outFiles[fn].write('%s\t%s\t%s\n' % (n,seq,qual))
inFile.close()
outFiles[0].close()
outFiles[1].close()

print 'Done'

