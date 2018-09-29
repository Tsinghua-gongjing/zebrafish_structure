#-*- coding:utf-8 -*-

from Environment import *

def calcScore(duplexGroup, supportReadFile, task_name, coverageMethod, genomeSizeFile=""):
    posBedFile = GLOBAL.outDir+"/tmp."+task_name+".cluster.pos.bed"
    negBedFile = GLOBAL.outDir+"/tmp."+task_name+".cluster.neg.bed"
    posGenomeCovFile = GLOBAL.outDir+"/tmp."+task_name+".pos.genomeCov"
    negGenomeCovFile = GLOBAL.outDir+"/tmp."+task_name+".neg.genomeCov"
    posIntCovFile = GLOBAL.outDir+"/tmp."+task_name+".pos.intCov"
    negIntCovFile = GLOBAL.outDir+"/tmp."+task_name+".neg.intCov"
    # define a Alert Object
    Alerter = ProcessAlert(total=len(duplexGroup), stepsToAlert=50000, OUT=GLOBAL.LOG, careRecentTimes=5)
    # write into pos/neg bed file
    dgCount = 0
    POS = open(posBedFile, 'w')
    NEG = open(negBedFile, 'w')
    for dg in sortObjDictByAttr(duplexGroup, ['chr1', 'strand1', 'chr2', 'strand2', 'start1', 'end1', 'start2', 'end2']):
        Alerter.alert()
        if duplexGroup[dg].collapsedTo != -1: continue
        dgCount += 1
        print >>(POS if duplexGroup[dg].strand1=='+' else NEG), "\t".join([duplexGroup[dg].chr1, `duplexGroup[dg].start1`, `duplexGroup[dg].end1`, `dg`, '1' , duplexGroup[dg].strand1])
        print >>(POS if duplexGroup[dg].strand2=='+' else NEG), "\t".join([duplexGroup[dg].chr2, `duplexGroup[dg].start2`, `duplexGroup[dg].end2`, `dg`, '2' , duplexGroup[dg].strand2])
    Alerter.finish()
    POS.close()
    NEG.close()
    if sortSupportParallel():
        #import multiprocessing
        #cores = multiprocessing.cpu_count()
        cmd_sort_sopport_read = "sort --parallel=8 -k1,1 -k2,3n %s | uniq > %s"
    else:
        cmd_sort_sopport_read = "sort -k1,1 -k2,3n %s | uniq > %s"
    runCMD(cmd_sort_sopport_read % (supportReadFile, supportReadFile+'.sorted'))
    if coverageMethod == 'pileup':
        #cmd_genomeCov = "awk 'NF==6{print $0}' %s | bedtools genomecov -i - -g %s -bg -strand %s > %s"
        cmd_genomeCov = "bedtools genomecov -i %s -g %s -bg -strand %s > %s"
        runCMD( cmd_genomeCov % (supportReadFile+'.sorted', genomeSizeFile, '+', posGenomeCovFile) )
        runCMD( cmd_genomeCov % (supportReadFile+'.sorted', genomeSizeFile, '-', negGenomeCovFile) )
        cmd_intersect = "bedtools intersect -wb -a %s -b %s > %s"
        runCMD( cmd_intersect % (posGenomeCovFile, posBedFile, posIntCovFile) )
        runCMD( cmd_intersect % (negGenomeCovFile, negBedFile, negIntCovFile) )
    elif coverageMethod == 'count':
        cmd_genomeCov = "bedtools intersect -a %s -b %s -c -s > %s"
        runCMD( cmd_genomeCov % (posBedFile, supportReadFile+'.sorted', posIntCovFile) )
        runCMD( cmd_genomeCov % (negBedFile, supportReadFile+'.sorted', negIntCovFile) )
    print >>GLOBAL.LOG, "Read file %s ..." % (posIntCovFile, )
    addCoverage ( posIntCovFile, duplexGroup, coverageMethod=coverageMethod );
    print >>GLOBAL.LOG, "Read file %s ..." % (negIntCovFile, )
    addCoverage ( negIntCovFile, duplexGroup, coverageMethod=coverageMethod );

def addCoverage(intCovFile, duplexGroup, coverageMethod):
    # define a Alert Object
    CovFileLines = countFileLines(intCovFile)
    Alerter = ProcessAlert(total=CovFileLines, stepsToAlert=50000, OUT=GLOBAL.LOG, careRecentTimes=10)
    # read
    ICF = open(intCovFile)
    line = ICF.readline()
    while line:
        Alerter.alert()
        data = line.strip().split()
        if coverageMethod == 'pileup':
            score = int(data[3])
            arm = int(data[8])
            dg = int(data[7])
        elif coverageMethod == 'count':
            score = int(data[6])
            arm = int(data[4])
            dg = int(data[3])
        if arm not in (1,2):
            print >>sys.stderr, "Unexpected Undefined Line: \n\t%s" % (line, )
            line = ICF.readline()
            continue
        if not defined(duplexGroup, dg):
            print >>sys.stderr, "Unexpected undefined duplexgroup: %i" % (dg, )
            line = ICF.readline()
            continue
        if duplexGroup[ dg ].cov == -1:
            duplexGroup[ dg ].cov = {1:0, 2:0}
            duplexGroup[ dg ].cov[ arm ] = score
        else:
            duplexGroup[ dg ].cov[ arm ] = max(score, duplexGroup[ dg ].cov[ arm ])
        line = ICF.readline()
    Alerter.finish()
    ICF.close()


# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

calcScore(duplexGroup, parameters['genomeSizeFile'], supportReadFile, task_name=parameters['task'], coverageMethod=parameters['coverage'])

# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 备份
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

import copy
duplexGroup_bak = copy.deepcopy(duplexGroup)

# 还原
duplexGroup = copy.deepcopy(duplexGroup_bak)

"""


