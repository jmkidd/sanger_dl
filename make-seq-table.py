import glob
import sys
from optparse import OptionParser

###############################################################################
def read_fasta_file_to_dict(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print 'ERROR, FILE DOESNNOT START WITH >'
        sys.exit()
    myName = line[1:].split()[0]
    myName = myName.split('|')[-1]
    
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:].split()[0]
            myName = myName.split('|')[-1]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
###############################################################################
def read_fasta_qual_file_to_dict(qualFile):
    myDict = {}
    inFile = open(qualFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print 'ERROR, FILE DOESNNOT START WITH >'
        sys.exit()
    myName = line[1:].split()[0]
    myName = myName.split('|')[-1]
    
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    myQualStr = ''
    

    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = myQualStr
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = myQualStr
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:].split()[0]
            myName = myName.split('|')[-1]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            myQualStr = ''
            continue
        qualList = line.split()        
        qs = [ chr(int(i) + 33) for i in qualList]
        qs = ''.join(qs)

        myQualStr += qs
        
    inFile.close()
    return myDict
###############################################################################
def int_to_qual(i):
    qoffSet = 33
    c = chr(i + qoffSet)
    return c
###############################################################################


USAGE = """
make-seq-table.py  --in <dir of input files>   

makes table of ti, seq, qual 1 per line for downstream processing

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inDir', help = 'input dir to process')


(options, args) = parser.parse_args()

if options.inDir is None:
    parser.error('inDir not given')

###############################################################################

if options.inDir[-1] != '/':
    options.inDir += '/'

faFiles = glob.glob(options.inDir + '*.fa')
print 'Found %i fa files' % len(faFiles)

outFileName = options.inDir + 'seq_qual.parsed'
print 'Writing output to',outFileName
outFile = open(outFileName,'w')
outFile.write('#ti\tseq\tqual\n')

for fa in faFiles:
    print fa
    qf = fa.replace('.fa','.qual')
    print qf
    faDict = read_fasta_file_to_dict(fa)
    print 'len fa',len(faDict)
    qaDict = read_fasta_qual_file_to_dict(qf)
    print 'len qa',len(qaDict)
        
    for name in faDict:
        outFile.write('%s\t%s\t%s\n' % (name,faDict[name]['seq'],qaDict[name]['seq']))
    

    
outFile.close()