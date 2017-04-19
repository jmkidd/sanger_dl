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







print 'Working on bam file',options.bam

