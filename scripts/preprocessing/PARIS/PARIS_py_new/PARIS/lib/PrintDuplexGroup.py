#-*- coding:utf-8 -*-

from Environment import *

def printDuplexGroup(outputFile, duplexGroup, read, readmap, minSupport=2, method="harmonic"):
    print >>GLOBAL.LOG, "Output Duplex Groups to file %s\n\t %s" % (outputFile, now())
    # define a Alert Object
    Alerter = ProcessAlert(total=len(duplexGroup), stepsToAlert=50000, OUT=GLOBAL.LOG, careRecentTimes=5)
    # output
    OUT = open(outputFile, 'w')
    dgCount = 0
    for dg in sortObjDictByAttr(duplexGroup, ['chr1', 'strand1', 'chr2', 'strand2', 'start1', 'end1', 'start2', 'end2']):
        Alerter.alert()
        if duplexGroup[dg].support < minSupport: continue
        if duplexGroup[dg].collapsedTo != -1: continue
        dgCount += 1
        left = max(duplexGroup[dg].cov[1], duplexGroup[dg].support) if duplexGroup[dg].cov != -1 else duplexGroup[dg].support
        right = max(duplexGroup[dg].cov[2], duplexGroup[dg].support) if duplexGroup[dg].cov != -1 else duplexGroup[dg].support
        score = supportScore ( duplexGroup[dg].support, left, right, method )
        OUT.writelines( "Group %d == position " % (dg, ) )
        OUT.writelines( "%s(%s):%d-%d|%s(%s):%d-%d" % (duplexGroup[dg].chr1, duplexGroup[dg].strand1, duplexGroup[dg].start1, duplexGroup[dg].end1, duplexGroup[dg].chr2, duplexGroup[dg].strand2, duplexGroup[dg].start2, duplexGroup[dg].end2) )
        OUT.writelines( ", support %d, left %d, right %d, score %f.\n---\n" % (duplexGroup[dg].support, left, right, score))
        reads = [ int(it) for it in  duplexGroup[dg].reads.split(';') ]
        for idx in range(len(reads)):
            names = read[ reads[idx] ].name.split(';')
            for name in names:
                print >>OUT, "\t%s\t%s|%s:%d-%d<=>%s|%s:%d-%d" % (name, read[ reads[idx] ].left.chr, read[ reads[idx] ].left.strand, read[ reads[idx] ].left.start, read[ reads[idx] ].left.end, read[ reads[idx] ].right.chr, read[ reads[idx] ].right.strand, read[ reads[idx] ].right.start, read[ reads[idx] ].right.end)
    Alerter.finish()
    OUT.close()
    print >>GLOBAL.LOG, "%d clusters output to file %s \n\tTime: %s" % (dgCount, outputFile, now())

def supportScore(support, left, right, method):
    import math
    score = -1
    if method == 'geometric':
        score = 1.0*(support + 1) / math.sqrt( (left+1)*(right+1) )
    elif method == 'harmonic':
        score = 1.0*(support + 1)*( 1.0/(left+1) + 1.0/(right+1) )/2
    return score


# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

printDuplexGroup(parameters['output'], duplexGroup, read, readmap, minSupport=parameters['minSupport'], method=parameters['scoringMethod'])

"""

