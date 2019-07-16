#!/usr/bin/env python
#-*- coding:utf-8 -*-

from lib.Environment import *
from lib.IntronInArea import *
from lib.GenPairClusterFromSamFile import *
from lib.GenPairClusterFromJunctionFile import *
from lib.GenDuplexGroup import *
from lib.CollapseDuplexGroup import *
from lib.FinalizeDuplexGroup import *
from lib.FinalizeReads import *
from lib.PrintSupportSam import *
from lib.CalcScore import *
from lib.PrintDuplexGroup import *


Usage = """

## --------------------------------------
 Call base pair groups from PARIS sequencing
 Command:
  %s -i input_sam_file -j chiastic_junction_file -s chiastic_support_sam_file -o output_read_group_fileName --genomeFile genome_fa --intronAnnoFile intronAnnotationFile

 # what it is:
  -i                    <String>
                                input sam file of spliced alignment, multiple files seperated by colon. example: Aligned1.out.sam:Aligned2.out.sam
  -j                    <String>
                                chiastic junction file, multiple files seperated by colon. example: Chimeric1.out.junction:Chimeric2.out.junction
  -s                    <String>
                                chiastic junction support alignment file, multiple files seperated by colon. example: Chimeric1.out.sam:Chimeric2.out.sam
  -o                    <String>
                                output file prefix (default: DG)
  --genomeFile          <String>
                                genome fasta file
  --intronAnnoFile      <String>
                                intron annotation file, no specify is OK

 # more options:
  * Temporary File:
  --tmpFileName         <String>
                                a tag to label all temporary files (defaule: random_number)

  * Logs and Errors:
  --log                 <String>
                                output log information to this file (default: stdout)
  --error               <String>
                                output error and warning information to this file (default: stderr)

  1. Generate Pair Cluster:
  --removeRedundancy    <yes/no>
                                remove redundancy in input sam or junctions (yes/no default: no)
  --minOverhang         <Interger>
                                minimum overhang length for valid mapping (default: 15)
  --localAlign          <yes/no>
                                local align both arms of a gapped read (yes/no default: no)
  --preserveMultimap    <yes/no>
                                preserve multi-mapped reads (yes/no default: no)
  --intronFlanking      <Interger>
                                allow how many nucleotides flack around the intron border (default: 3)

  2. Generate Duplex Group:
  --minOverlap          <Interger>
                                minimum number of overlap nucleotides to define a read duplex (default: 5)
  --multipleDG          <yes/no>
                                allow one read in multiple duplex groups (yes/no default: no)

  3. Collapse Duplex Group:
  --maxGap              <Interger>
                                how many nucleotides allowed gapped during merging duplex group (default: 10)
  --maxDGOverhang       <Interger>
                                how long the DG arm allowed (default: 30)

  4. Scoring and Filter Duplex Group:
  --coverage            <pileup/count>
                                estimate read count in a duplex arm by pileup or count. if specify pileup, --genomeSizeFile is neccessary (default: pileup)
  --genomeSizeFile      <String>        
                                genome size file
  --minSupport          <Interger>      
                                minimum number of supporting reads (default: 2)
  --scoringMethod       <harmonic/geometric>
                                scoring method (default: harmonic)

""" % (sys.argv[0], )

def init():
    GLOBAL.init_GLOBAL()
    params = {
        "output": "DG", "tmpFileName": str(random.randint(10000,30000)), "removeRedundancy": False, "minOverhang": 15, \
        "localAlign": False, "preserveMultimap": False, "intronFlanking": 3, "minOverlap": 5, "multipleDG":False, \
        "maxGap": 10, "maxDGOverhang": 30, "coverage": "pileup", "minSupport": 2, "scoringMethod": "harmonic"
    }
    opts, args = getopt.getopt(sys.argv[1:], 'hi:j:s:o:', \
        ['genomeFile=', 'intronAnnoFile=', 'tmpFileName=', 'log=', 'error=', 'removeRedundancy=', 'minOverhang=', 'localAlign=', 'preserveMultimap=', 'intronFlanking=', 'minOverlap=', 'multipleDG=', 'maxGap=', 'maxDGOverhang=','coverage=', 'genomeSizeFile=', 'minSupport=', 'scoringMethod='])
    for op, value in opts:
        if op == '-h':
            print >>sys.stdout, Usage;
            sys.exit(-1)
        # Basic Parameters
        elif op == '-i':
            params['samFiles'] = os.path.abspath(value)
            #assert os.path.exists(params['samFiles']), "Error: -i file doesn't exist"
        elif op == '-j':
            params['chiasticFiles'] = os.path.abspath(value)
            #assert os.path.exists(params['chiasticFiles']), "Error: -j file doesn't exist"
        elif op == '-s':
            params['chiasticSamFiles'] = os.path.abspath(value)
            #assert os.path.exists(params['chiasticSamFiles']), "Error: -s file doesn't exist"
        elif op == '-o':
            params['output'] = os.path.abspath(value)
            if '/' in params['output']:
                assert os.path.isdir( os.path.dirname(params['output']) ), "Error: -o is not a directory"
            GLOBAL.outDir = os.path.dirname(params['output'])
        elif op == '--genomeFile':
            params['genomeFile'] = os.path.abspath(value)
            assert os.path.exists(params['genomeFile']), "Error: --genomeFile file doesn't exist"
        elif op == '--intronAnnoFile':
            params['intronAnnoFile'] = os.path.abspath(value)
            assert os.path.exists(params['intronAnnoFile']), "Error: --intronAnnoFile file doesn't exist"
        # Temporary File
        elif op == '--tmpFileName':
            params['tmpFileName'] = value
        # Logs and Errors
        elif op == '--log':
            params['log'] = os.path.abspath(value)
            GLOBAL.LOG = OpenFile( params['log'], 'w')
        elif op == '--error':
            params['error'] = os.path.abspath(value)
            GLOBAL.ERR = OpenFile( params['error'], 'w' )
        # 1. Generate Pair Cluster
        elif op == '--removeRedundancy':
            assert value in ('yes', 'no'), "Error: --removeRedundancy must be yes or no"
            params['removeRedundancy'] = True if value == "yes" else False
        elif op == '--minOverhang':
            params['minOverhang'] = int(value)
        elif op == '--localAlign':
            assert value in ('yes', 'no'), "Error: --localAlign must be yes or no"
            params['useLongestPairing'] = True if value == "yes" else False
        elif op == '--preserveMultimap':
            assert value in ('yes', 'no'), "Error: --preserveMultimap must be yes or no"
            params['preserveMultimap'] = True if value == "yes" else False
        elif op == '--intronFlanking':
            params['intronFlanking'] = int(value)
        # 2. Generate Duplex Group
        elif op == '--minOverlap':
            params['minOverlap'] = int(value)
        elif op == '--multipleDG':
            assert value in ('yes', 'no'), "Error: --multipleDG must be yes or no"
            params['multipleDG'] = True if value == "yes" else False
        # 3. Collapse Duplex Group
        elif op == '--maxGap':
            params['maxGap'] = int(value)
        elif op == '--maxDGOverhang':
            params['maxDGOverhang'] = int(value)
        # 4. Scoring and Filter Duplex Group
        elif op == '--coverage':
            params['coverage'] = value
            assert params['coverage'] in ('pileup', 'count'), "Error: --coverage must be one of pileup/count"
        elif op == '--genomeSizeFile':
            params['genomeSizeFile'] = os.path.abspath(value)
            assert os.path.exists(params['genomeSizeFile']), "Error: --genomeSizeFile file doesn't exist"
        elif op == '--minSupport':
            params['minSupport'] = int(value)
        elif op == '--scoringMethod':
            params['scoringMethod'] = value
            assert params['scoringMethod'] in ('harmonic', 'geometric'), "Error: --scoringMethod must be one of harmonic/geometric"
    if ('samFiles' not in params) or ('samFiles' not in params and 'chiasticSamFiles' not in params) or ('genomeFile' not in params):
          print >>sys.stdout, Usage
          sys.exit(-1)
    if (params['coverage'] == "pileup") and ('genomeSizeFile' not in params):
        print >>sys.stdout, "Error: pipeup coverage method need a <<genome size file>>\n\t%s" % (now(), )
        print >>sys.stdout, Usage
        sys.exit(-1)
    return params


def beautifulPrintParams(Parameters):
    print >>GLOBAL.LOG, "\n# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
    print >>GLOBAL.LOG, "                    Show Parameters                        "
    for Key in Parameters.keys():
        print >>GLOBAL.LOG, "%25s\t%s" % (Key, Parameters[Key])
    print >>GLOBAL.LOG, "# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"
    flush()

StepTags = ['1. genPairClusterFromSamFile', '2. genPairClusterFromJunctionFile', '3. genDuplexGroup', '4. collapseDuplexGroup', '5. finalizeDuplexGroup', '6. finalizeReads', '7. printSupportSam', '8. calcScore', '9. printDuplexGroup']

def main():
    parameters = init()
    beautifulPrintParams(parameters)
    if 'genomeFile' in parameters:
        showTime('read genome file: %s' % (parameters['genomeFile'], ))
        #global Genome
        GLOBAL.Genome = SeqFunc.seqClass(parameters['genomeFile'])
    if 'intronAnnoFile' in parameters:
        showTime('read intron: %s' % (parameters['intronAnnoFile'], ))
        #global Intron
        GLOBAL.Intron = intronInAreaClass(parameters['intronAnnoFile'])
    # define global variable
    read = {}; readmap = {}; duplexGroup = {}
    allSupportSam = ''
    if 'samFiles' in parameters:
        #           ==================1. genPairClusterFromSamFile==================
        print >>GLOBAL.LOG, "="*20+StepTags[0]+"="*20
        print >>GLOBAL.ERR, "="*20+StepTags[0]+"="*20
        allSupportSam = parameters['samFiles']
        samFiles = parameters['samFiles'].split(':')
        for samFile in samFiles:
            showTime('genPairClusterFromSamFile: %s' % (samFile, ))
            genPairClusterFromSamFile( samFile, read, readmap, parameters['tmpFileName'], \
                removeRedundancy=parameters['removeRedundancy'], \
                minOverhang=parameters['minOverhang'], \
                useLongestPairing=parameters['useLongestPairing'], \
                preserveMultimap=parameters['preserveMultimap'], \
                intronFlanking=parameters['intronFlanking'], \
                verbose=False)
    if 'chiasticFiles' in parameters:
        #           ==================2. genPairClusterFromJunctionFile==================
        print >>GLOBAL.LOG, "="*20+StepTags[1]+"="*20
        print >>GLOBAL.ERR, "="*20+StepTags[1]+"="*20
        if allSupportSam:
            allSupportSam += ':' + parameters['chiasticSamFiles']
        else:
            allSupportSam += parameters['chiasticSamFiles']
        chiasticFiles = parameters['chiasticFiles'].split(':')
        for chiasticFile in chiasticFiles:
            showTime('genPairClusterFromJunctionFile: %s' % (chiasticFile, ))
            genPairClusterFromJunctionFile ( chiasticFile, read, readmap, parameters['tmpFileName'], \
                removeRedundancy=parameters['removeRedundancy'], \
                minOverhang=parameters['minOverhang'],\
                useLongestPairing=parameters['useLongestPairing'], \
                intronFlanking=parameters['intronFlanking'], \
                verbose=False )
    #sys.exit(-1)
    #       ==================3. genDuplexGroup==================
    print >>GLOBAL.LOG, "="*20+StepTags[2]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[2]+"="*20
    genDuplexGroup(duplexGroup, read, task_name=parameters['tmpFileName'], minOverlap=parameters['minOverlap'], multipleDG=parameters['multipleDG'])
    #       ==================4. collapseDuplexGroup==================
    print >>GLOBAL.LOG, "="*20+StepTags[3]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[3]+"="*20
    collapseDuplexGroup(duplexGroup, read, maxGap=parameters['maxGap'], maxTotal=parameters['maxDGOverhang'])
    #       ==================5. finalizeDuplexGroup==================
    print >>GLOBAL.LOG, "="*20+StepTags[4]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[4]+"="*20
    finalizeDuplexGroup(duplexGroup, read, readmap, task_name=parameters['tmpFileName'])
    #       ==================6. finalizeReads==================
    print >>GLOBAL.LOG, "="*20+StepTags[5]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[5]+"="*20
    old_finalizeReads(read, readmap, duplexGroup)
    #       ==================7. printSupportSam==================
    print >>GLOBAL.LOG, "="*20+StepTags[6]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[6]+"="*20
    supportSamFile = parameters["output"]+'.sam'; supportReadFile = parameters["output"]+'.reads'
    printSupportSam(supportSamFile, allSupportSam, read, readmap, outputReadFile=supportReadFile)
    #       ==================8. calcScore==================
    print >>GLOBAL.LOG, "="*20+StepTags[7]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[7]+"="*20
    calcScore(duplexGroup, supportReadFile, task_name=parameters['tmpFileName'], coverageMethod=parameters['coverage'], genomeSizeFile=parameters['genomeSizeFile'])
    #       ==================9. calcScore==================
    print >>GLOBAL.LOG, "="*20+StepTags[8]+"="*20
    print >>GLOBAL.ERR, "="*20+StepTags[8]+"="*20
    printDuplexGroup(parameters['output'], duplexGroup, read, readmap, minSupport=parameters['minSupport'], method=parameters['scoringMethod'])


main()














