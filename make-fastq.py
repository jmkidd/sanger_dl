import glob
import sys
import gzip

from optparse import OptionParser
###############################################################################
USAGE = """
make-fastq.py  --in <dir>   

will make fq.gz, run after doing the sort

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inDir', help = 'input dir to process')


(options, args) = parser.parse_args()

if options.inDir is None:
    parser.error('inDir not given')

###############################################################################

if options.inDir[-1] != '/':
    options.inDir += '/'



baseName = options.inDir.split('/')[-2]

fq1Name= options.inDir + baseName + '.R1.fq.gz'
fq1Table = options.inDir + 'seq_qual.R1.sorted'

print 'Reading %s writing %s' % (fq1Table,fq1Name)
inFile = open(fq1Table,'r')
outFile = gzip.open(fq1Name,'w')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    outFile.write('@%s\n%s\n+\n%s\n' % (line[0],line[1],line[2]))

outFile.close()
inFile.close()


fq2Name= options.inDir + baseName + '.R2.fq.gz'
fq2Table = options.inDir + 'seq_qual.R2.sorted'

print 'Reading %s writing %s' % (fq2Table,fq2Name)
inFile = open(fq2Table,'r')
outFile = gzip.open(fq2Name,'w')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    outFile.write('@%s\n%s\n+\n%s\n' % (line[0],line[1],line[2]))

outFile.close()
inFile.close()