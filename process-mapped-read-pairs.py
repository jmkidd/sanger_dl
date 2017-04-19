import sys
import os
import genutils
import pysam

from optparse import  OptionParser
#############################################################################

USAGE = """process-mapped-read-pairs.py  --bam <bam file> 

input must be a query-name-sorted BAM file.

Will process and print out various things.




"""
parser = OptionParser(USAGE)
parser.add_option('--bam',dest='bam', help = 'bam name, query sorted')

(options,args)=parser.parse_args()

if options.bam is None:
    parser.error('bam  not given')
#############################################################################
def process_read_pair(r1,r2,outFile):
#    print r1.query_name
#    print r2.query_name
    # check for mark dups
    if r1.is_duplicate is True:
        return False
    if r2.is_duplicate is True:
        return False
        
        
    print r1
    print r2
    sys.exit()


    return True


#############################################################################




print 'Working on bam file',options.bam

samFile = pysam.AlignmentFile(options.bam,'r')
outFileName = options.bam + '.readpairtable'
outFile = open(outFileName,'w')

readCache = [] # cache where we will keep the read records for processing

numPairs = 0
numPairsPass = 0
for read in samFile:
    #skip over not primary alignments
    if read.is_secondary is True:
        continue
    if read.is_supplementary is True:
        continue
    
    if len(readCache) == 0:
        readCache.append(read)
    elif len(readCache) == 1:
        if read.query_name == readCache[0].query_name:
            numPairs += 1
            res = process_read_pair(readCache[0],read,outFile)
            if res is True:
                numPairsPass += 1
            readCache.pop()
        else:
            print 'CACHE not match!!'
            print readCache[0]
            print read
            print 'CACHE not match!!'
            sys.exit()
    else:
        print 'len readcache is',len(readCache)
        print readCache
        print read
        print 'len readcache is',len(readCache)
        sys.exit()

samFile.close()    
print 'Found %i read pairs' % numPairs
print 'Found %i read pairs pass ' % numPairsPass

outFile.close()