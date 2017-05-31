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
    t = ['.','.','.','.','.','.','.','.','.']


#    print r1.query_name
#    print r2.query_name
    # check for mark dups
    if r1.is_duplicate is True:
        return False
    if r2.is_duplicate is True:
        return False
        
        
#    print r1
    # at this point, we have an alignment to consider
    # build a list of alignments of structure [qPos,rPos,qBp,rBp]
    # with '-' bp for gaps in query/reference.  length comes from
    # get_aligned_pairs
    
    if r1.is_unmapped is False:
        readSeq = r1.query_sequence
        quals = r1.query_qualities    
        alignedSeq = r1.get_aligned_pairs(with_seq=True)
        alignData = []
        # setup empty list of 4-lits.  could expand this to consider read qual
        for i in range(len(alignedSeq)):
            alignData.append(['.','.','.','.'])
        # go through and populate with appropriate info
        for i in range(len(alignedSeq)):
            alignData[i][0] = alignedSeq[i][0] 
            alignData[i][1] = alignedSeq[i][1]
            if alignData[i][0] == None:
                 alignData[i][2] = '-' 
            else:
                 alignData[i][2] = readSeq[alignData[i][0]].upper() #pull out read Bp from readseq
    
            if alignData[i][1] == None:
                alignData[i][3] = '-'
            else:
                alignData[i][3] = alignedSeq[i][2].upper()
        # alignData structure is now built and ready to parse!
        R1mapRes = calc_match_mismatch(alignData,quals)
        refPos = r1.get_reference_positions()
        r1Start = refPos[0] + 1
        r1End = refPos[-1] + 1

    if r2.is_unmapped is False:
        readSeq = r2.query_sequence
        quals = r2.query_qualities    
        alignedSeq = r2.get_aligned_pairs(with_seq=True)
        alignData = []
        # setup empty list of 4-lits.  could expand this to consider read qual
        for i in range(len(alignedSeq)):
            alignData.append(['.','.','.','.'])
        # go through and populate with appropriate info
        for i in range(len(alignedSeq)):
            alignData[i][0] = alignedSeq[i][0] 
            alignData[i][1] = alignedSeq[i][1]
            if alignData[i][0] == None:
                 alignData[i][2] = '-' 
            else:
                 alignData[i][2] = readSeq[alignData[i][0]].upper() #pull out read Bp from readseq
    
            if alignData[i][1] == None:
                alignData[i][3] = '-'
            else:
                alignData[i][3] = alignedSeq[i][2].upper()
        # alignData structure is now built and ready to parse!
        R2mapRes = calc_match_mismatch(alignData,quals)
        refPos = r2.get_reference_positions()
        r2Start = refPos[0] + 1
        r2End = refPos[-1] + 1
    #build our list of clone assignments
    
    nl = []
    nl.append(r1.query_name)
    if r1.is_unmapped is False and r2.is_unmapped is False: # both are mapped
        if r1Start <= r2Start:  # r1 goes first
            nl.append(r1.reference_name)
            nl.append(r1Start)
            nl.append(r1End)
            nl.append(r1.mapping_quality)
            if r1.is_reverse is True:
                nl.append('-')
            else:
                nl.append('+')
            nl.extend(R1mapRes)            
            #now r2
            nl.append(r2.reference_name)
            nl.append(r2Start)
            nl.append(r2End)
            nl.append(r2.mapping_quality)
            if r2.is_reverse is True:
                nl.append('-')
            else:
                nl.append('+')
            nl.extend(R2mapRes)
        else:  #r2 goes first
            nl.append(r2.reference_name)
            nl.append(r2Start)
            nl.append(r2End)
            nl.append(r2.mapping_quality)
            if r2.is_reverse is True:
                nl.append('-')
            else:
                nl.append('+')
            nl.extend(R1mapRes)            
            #now r1
            nl.append(r1.reference_name)
            nl.append(r1Start)
            nl.append(r1End)
            nl.append(r1.mapping_quality)
            if r1.is_reverse is True:
                nl.append('-')
            else:
                nl.append('+')
            nl.extend(R1mapRes)
    elif r1.is_unmapped is True and r2.is_unmapped is True: # both unmapped
        nl.extend(t)
        nl.extend(t)
    elif  r1.is_unmapped is False: # r1 mapped, r2 not
        nl.append(r1.reference_name)
        nl.append(r1Start)
        nl.append(r1End)
        nl.append(r1.mapping_quality)
        if r1.is_reverse is True:
            nl.append('-')
        else:
            nl.append('+')
        nl.extend(R1mapRes)            
        nl.extend(t)
    else: #r2 mapped, r1 not                
        nl.append(r2.reference_name)
        nl.append(r2Start)
        nl.append(r2End)
        nl.append(r2.mapping_quality)
        if r2.is_reverse is True:
            nl.append('-')
        else:
            nl.append('+')
        nl.extend(R2mapRes)            
        nl.extend(t)
    # ready to print
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    return True
#############################################################################
def calc_match_mismatch(alignData,quals):
    numMatch = 0
    numMismatch = 0
    numMatchQ30 = 0
    numMistmatchQ30 = 0
    for i in range(len(quals)):
        if alignData[i][2] == '-' or alignData[i][3] == '-':
            continue         
        if alignData[i][2] == alignData[i][3]:
            numMatch += 1
            if quals[alignData[i][0]] >= 30:
                numMatchQ30 += 1
        else:
            numMismatch += 1
            if quals[alignData[i][0]] >= 30:
                numMistmatchQ30 += 1
    tot = numMatch + numMismatch
    totQ30 = numMatchQ30 + numMistmatchQ30
    if tot >0:
        tf = float(numMatch)/tot
    else:
        tf = 0.0

    if totQ30 >0:
        qtf = float(numMatchQ30)/totQ30
    else:
        qtf = 0.0
    
    return [tot,tf,totQ30,qtf]
#############################################################################

print 'Working on bam file',options.bam

samFile = pysam.AlignmentFile(options.bam,'r')
outFileName = options.bam + '.readpairtable'
outFile = open(outFileName,'w')

nl = ['#cloneName','chrom','start','end','mapQ','dir','alignedBp','fracMatch','alignedBpQ30','fracMatchQ30','chrom','start','end','dir','alignedBp','fracMatch','alignedBpQ30','fracMatchQ30']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

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
    if numPairs % 5000 == 0:
        print 'Did %i pairs...' % numPairs
samFile.close()    
print 'Found %i read pairs' % numPairs
print 'Found %i read pairs pass ' % numPairsPass
outFile.close()