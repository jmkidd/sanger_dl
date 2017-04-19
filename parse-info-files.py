import glob
import sys
from optparse import OptionParser


###############################################################################
#parse info file into dictionary
def get_info_record(inFile):
    rec = {}
    while True:
        line = inFile.readline()
        if line == '': #end of file?
           return {} 
        if line == '\n': #end of record
            return rec
        line = line.rstrip()
        line = line.split('\t')
        n = line[0]
        v = line[1]
        n = n.replace(':','')
        rec[n] = v
###############################################################################
###############################################################################
USAGE = """
parse-info-files.py  --in <dir of input files>   

will parse info files to get name, info, etc from runs.
"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inDir', help = 'input dir to process')


(options, args) = parser.parse_args()

if options.inDir is None:
    parser.error('inDir not given')

###############################################################################

if options.inDir[-1] != '/':
    options.inDir += '/'

nameType = 2 # type 1: just take clone ID -- this is right info, for example some BACs
             # type 2: combine plate and well to clone ID

infoFiles = glob.glob(options.inDir + '*.info')
print 'Found %i info files' % len(infoFiles)
outFileName = options.inDir + 'info.parsed'
print 'Writing output to',outFileName
outFile = open(outFileName,'w')
outFile.write('#ti\tclone_wellID\ttrace_end\n')
for fn in infoFiles:
    print fn
    inFile = open(fn,'r')
    while True:
        rec = get_info_record(inFile)
        if rec == {}:
            break
        if nameType == 1:
            cloneWell = rec['clone_id']
        if nameType == 2:
            cloneWell = rec['plate_id'] + '_' + rec['well_id']

        outFile.write('%s\t%s\t%s\n' % (rec['ti'],cloneWell,rec['trace_end']))        
    inFile.close()
    
outFile.close()