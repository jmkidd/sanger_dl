import glob
import sys
from optparse import OptionParser
###############################################################################
USAGE = """
count-paired-from-info.py  --in <parsed info file>   

will count number of records, and number that are paired properly
"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFile', help = 'input dir to process')

(options, args) = parser.parse_args()

if options.inFile is None:
    parser.error('inFile not given')

###############################################################################
###############################################################################

numReads = 0

n2i = {} # name to index
n2i['REVERSE'] = 1
n2i['FORWARD'] = 0


inFile = open(options.inFile,'r')
cloneRecords = {}
dirCounts = [0,0]
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split('\t')
    numReads += 1
    n = line[1]
    d = line[2]
    i = n2i[d]
    dirCounts[i] += 1
    if n not in cloneRecords:
        cloneRecords[n] = [0,0]
    cloneRecords[n][i] += 1
inFile.close()
print '\nFile',options.inFile
print 'Number of reads:',numReads
print 'Forward %i reverse %i' % (dirCounts[0],dirCounts[1])
print 'Total clones',len(cloneRecords)

numWith1Fand1R = 0
# get clones with 1 forward and 1 reverse

outFilePairedFile = options.inFile + '.paired'
outFile = open(outFilePairedFile,'w')
for n in cloneRecords:
    if cloneRecords[n][0] == 1 and cloneRecords[n][1] == 1:
        numWith1Fand1R += 1
        outFile.write('%s\n' % (n))
outFile.close()
print 'Clones with 1 forward and 1 reverse',numWith1Fand1R
print 'Written to',outFilePairedFile