import requests
import sys
import math
import os
from datetime import datetime

from optparse import OptionParser


###############################################################################
USAGE = """
download-trace-files.py  --in <query>   --out <output dir>

will download fasta, qual, and info files into output dir as specificed

If all three files exist for a given 40,000 "page" that page will be skipped
This permits restart of interrupted downloads.

Count reocrds at the end to be sure all expected data was obtained.

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='query', help = 'input query to download')
parser.add_option('--out',dest='outDir', help = 'dir for output')

(options, args) = parser.parse_args()

if options.query is None:
    parser.error('query not given')
if options.outDir is None:
    parser.error('outDir not given')

###############################################################################
if options.outDir[-1] != '/':
    options.outDir += '/'

url = 'https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=raw'
recordsPerSet = 40000
#first lets figure out how many elements there are


q = 'query count ' + options.query
print 'Doing query to find number of elements'
print q
print str(datetime.now())
payload = {'query':q}
r = requests.post(url, data=payload)
numRecords = int(r.text)
print 'There are %i records to download' % numRecords

numSets = int(math.ceil(float(numRecords)/float(recordsPerSet)))
print 'How many sets needed?', float(numRecords)/float(recordsPerSet),numSets

for i in range(numSets):
    print '%i of %i' % (i,numSets-1)
    recNamesFile = options.outDir + '%s.names.txt' % i
    recFastaFile = options.outDir + '%s.fa' % i
    recQualFile = options.outDir + '%s.qual' % i
    recInfoFile = options.outDir + '%s.info' % i

    # check to see if all exists, if so do not re-download
    # this is useful for when redoing interrupted downloads
    if os.path.isfile(recNamesFile) is True and \
    os.path.isfile(recFastaFile) is True and \
    os.path.isfile(recQualFile) is True and \
    os.path.isfile(recInfoFile) is True: 
        continue
    


    print str(datetime.now())    
    
    #get the set IDs
    q = 'query page_size %i page_number %i %s' % (recordsPerSet,i,options.query)
    print q
    payload = {'query':q}
    r = requests.post(url, data=payload)
    outFile = open(recNamesFile,'w')
    outFile.write(r.text)
    outFile.close()
    tiList = []
    inFile = open(recNamesFile,'r')
    for ti in inFile:
        ti = ti.rstrip()
        ti = ti.lstrip()
        tiList.append(ti)
    inFile.close()
    print 'Have %i TI' % len(tiList)
    tiStr = ','.join(tiList)
    
    # get fasta
    q = 'retrieve fasta %s' % tiStr
    print 'Downloading fasta ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    outFile = open(recFastaFile,'w')
    outFile.write(r.text)    
    outFile.close()

    # get qual
    q = 'retrieve quality %s' % tiStr
    print 'Downloading quality ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    outFile = open(recQualFile,'w')
    outFile.write(r.text)    
    outFile.close()

    # get info
    q = 'retrieve info %s' % tiStr
    print 'Downloading info ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    outFile = open(recInfoFile,'w')
    outFile.write(r.text)    
    outFile.close()

print 'All done!'
print str(datetime.now())    
    

    





    
    

    
    




