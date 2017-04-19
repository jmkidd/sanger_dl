import requests
import sys
import math
from datetime import datetime


tiList = []
inFile = open('CH82.endinfo_9615.out','r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split('\t')
    ti = line[3]
    tiList.append(ti)
inFile.close()

print 'Read in %i TI' % len(tiList)

fNum = 0
for i in range(0,len(tiList),40000):
    print i,i+40000
    outFile = open('CHORI-82/%i.names.txt' % fNum,'w')
    nl = tiList[i:i+40000]
    nl = '\n'.join(nl) + '\n'
    outFile.write(nl)    
    outFile.close()
    fNum += 1
    
print 'fNum is',fNum
url = 'https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=raw'
recordsPerSet = 40000


for i in range(fNum):
    print '%i of %i' % (i,fNum-1)
    print str(datetime.now())    
    tiFile =  open('CHORI-82/%i.names.txt' % i,'r')
    tiList = []
    for line in tiFile:
        tiList.append(line)
    tiFile.close()
    tiStr = ','.join(tiList)
    
    # get fasta
    q = 'retrieve fasta %s' % tiStr
    print 'Downloading fasta ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    recFastaFile = 'CHORI-82/' + '%s.fa' % i
    outFile = open(recFastaFile,'w')
    outFile.write(r.text)    
    outFile.close()

    # get qual
    q = 'retrieve quality %s' % tiStr
    print 'Downloading quality ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    recQualFile = 'CHORI-82/' + '%s.qual' % i
    outFile = open(recQualFile,'w')
    outFile.write(r.text)    
    outFile.close()

    # get info
    q = 'retrieve info %s' % tiStr
    print 'Downloading info ...'
    print str(datetime.now())    

    payload = {'query':q}
    r = requests.post(url, data=payload)
    recInfoFile = 'CHORI-82/' + '%s.info' % i
    outFile = open(recInfoFile,'w')
    outFile.write(r.text)    
    outFile.close()

print 'All done!'
print str(datetime.now())    

    
    