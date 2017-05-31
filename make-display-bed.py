import sys
import os
import numpy as np

from optparse import  OptionParser
#############################################################################

USAGE = """make-display-bed.py  --in <CHORI-82.sel.byname.bam.readpairtable> 

input must be read-pair table as made by process-mapped-read-pairs.py

Will process and print out various things.


"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFile', help = 'file table name')

(options,args)=parser.parse_args()

if options.inFile is None:
    parser.error('in file  not given')
#############################################################################

print 'infile is',options.inFile

cloneListInfo = []
inFile = open(options.inFile,'r')

pairedList = []
inverList = []
eversionList = []
splitChromList = []

for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split('\t')

    cloneName = line[0]
    c1 = line[1]
    if c1 != '.':
       s1 = int(line[2])
       e1 = int(line[3])
       dir1 = line[5]
    else:
       c1 = '.'
       s1 = '.'
       e1 = '.'
       dir1 = '.'

 
    c2 = line[10]
    if c2 != '.':
        s2 = int(line[11])
        e2 = int(line[12])
        dir2 = line[14]
    else:
        s2 = '.'
        e2 = '.'
        dir2 = '.'
    
    # hacky deal with chrUn for now...
    if c1 == 'chrUn' or c2 == 'chrUn':
        continue
    
    
    #proper pairs:
    if c1 == '.' and c2 == '.':
        continue
    if c1 == c2 and dir1 == '+' and dir2 == '-':
        pairedList.append([cloneName,c1,s1,e1,s2,e2,'+-'])
    elif c1 == c2 and dir1 == '+' and dir2 == '+':
        inverList.append([cloneName,c1,s1,e1,s2,e2,'++'])
    elif c1 == c2 and dir1 == '-' and dir2 == '-':
        inverList.append([cloneName,c1,s1,e1,s2,e2,'--'])
    elif c1 == c2 and dir1 == '-' and dir2 == '+':
        eversionList.append([cloneName,c1,s1,e1,s2,e2,'-+'])
    elif c1 != c2:
        if c1 != '.':
            splitChromList.append([cloneName,c1,s1,e1,dir1])
        if c2 != '.':
            splitChromList.append([cloneName,c2,s2,e2,dir2])        
    
    else:
         print line
         sys.exit()
    
inFile.close()

print 'Read in !'
print 'Paired list',len(pairedList)
print 'inverList',len(inverList)
print 'eversionList',len(eversionList)
print 'splitChrom',len(splitChromList)

print 'Going throught paired list to get sizes!'
sizeList = []
for i in pairedList:
    s = i[5] - i[2] + 1
    sizeList.append(s)

print len(sizeList)

meanSize = np.mean(sizeList)
medianSize = np.median(sizeList)
stdSize = np.std(sizeList)

mad = np.median(np.abs(sizeList - medianSize))

print 'Mean %.1f median %.1s std %.1f' % (meanSize,medianSize,stdSize)
print 'MAD %.1f' % mad
lc = medianSize - 5*mad
uc = medianSize + 5*mad
print '5*mad cutoffs: %.1f - %.1f' % (lc,uc)

numInSize = 0
for i in sizeList:
    if i >= lc and i <= uc:
        numInSize += 1
print 'In size range: %i  %.3f of those considered' % (numInSize, float(numInSize)/len(sizeList))


# print out the concordant BAC clones...
concordantFileName = options.inFile + '.concordant.bed'
outFile = open(concordantFileName,'w')
for s in pairedList:
    sz = s[5] - s[2] + 1
    if sz >= lc and sz <= uc:
        nl = [s[1],s[2]-1,s[5],s[0],1000,'+',s[2]-1,s[5],'0,0,0',2]
        sz1 = s[3]-s[2] + 1
        sz2 = s[5]-s[4] + 1
        nl.append('%i,%i,' % (sz1,sz2))
        st1 = s[2] - s[2]
        st2 = s[4] - s[2]
        nl.append('%i,%i' % (st1,st2))
        nl = [str(j) for j in nl]
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)
        
outFile.close()


